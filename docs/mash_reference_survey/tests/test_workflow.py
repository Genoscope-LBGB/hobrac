"""Real local Snakemake DAG tests with fake hobrac/Mash and a loopback FASTA server."""

import gzip
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import threading
import time
import unittest
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from ..scripts import survey
from . import test_survey


@unittest.skipUnless(
    shutil.which("snakemake") and importlib.util.find_spec("matplotlib"),
    "Install Snakemake >=9, its Slurm executor plugin, and matplotlib",
)
class WorkflowTests(unittest.TestCase):
    def setUp(self):
        fixture = test_survey.SurveyTests()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        self.work = fixture.work
        self.rows = fixture.prepare()
        self.events = self.work / "events.jsonl"
        events = self.events
        counter = []

        class Handler(BaseHTTPRequestHandler):
            def do_GET(self):
                accession = self.path.strip("/")
                counter.append(accession)
                # One fast download lets Mash start while another is still downloading.
                delay = 0.05 if len(counter) == 1 else 1.5
                with events.open("a") as handle:
                    handle.write(
                        json.dumps(
                            {
                                "kind": "download",
                                "event": "start",
                                "accession": accession,
                                "time": time.monotonic(),
                            }
                        )
                        + "\n"
                    )
                time.sleep(delay)
                data = gzip.compress(b">chr1\nACGTACGTACGT\n")
                self.send_response(200)
                self.send_header("Content-Length", str(len(data)))
                self.end_headers()
                self.wfile.write(data)
                with events.open("a") as handle:
                    handle.write(
                        json.dumps(
                            {
                                "kind": "download",
                                "event": "end",
                                "accession": accession,
                                "time": time.monotonic(),
                            }
                        )
                        + "\n"
                    )

            def log_message(self, format, *args):
                pass

        server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        self.addCleanup(server.server_close)
        self.addCleanup(server.shutdown)
        for row in self.rows:
            row["fasta_url"] = (
                f"http://127.0.0.1:{server.server_port}/{row['accession']}"
            )
        survey.write_tsv(
            self.work / "manifest.tsv", self.rows, ("task_id", *survey.FIELDS)
        )
        snapshot = json.loads((self.work / "snapshot.json").read_text())
        snapshot["files"]["manifest.tsv"] = survey.sha256(self.work / "manifest.tsv")
        survey.save_json(self.work / "snapshot.json", snapshot)
        self.scripts = Path(__file__).resolve().parents[1]
        self.config = self.work / "config.sh"
        self.config.write_text(
            (self.scripts / "config.example.sh")
            .read_text()
            .replace("/path/to/shared/scratch/hobrac_mash_survey", str(self.work))
            .replace("MAX_CONCURRENT=4", "MAX_CONCURRENT=2")
            .replace("CPUS=2", "CPUS=1")
        )
        bin_dir = self.work / "bin"
        bin_dir.mkdir()
        hobrac = bin_dir / "hobrac"
        hobrac.write_text(
            f"#!{sys.executable}\n"
            + """
import json, os, sys, time
from pathlib import Path
args = sys.argv[1:]
assert args[args.index("-e") + 1] == "local"
assert "--stop-after-mash" in args
assert "--rerun-incomplete" in args
assert Path(os.environ["TAXONKIT_DB"]).is_dir()
output = Path(args[args.index("-o") + 1])
accession = output.parent.name
assembly = output / "assembly"
assembly.mkdir(exist_ok=True)
(assembly / (accession + ".fna")).write_text(">chr1\\nACGT\\n")
(assembly / (accession + ".chr_rename.tsv")).write_text("old_name\\tnew_name\\n")
def event(phase):
    with open(os.environ["TEST_EVENTS"], "a") as handle:
        record = {"kind": "mash", "event": phase,
                  "accession": accession, "time": time.monotonic()}
        handle.write(json.dumps(record) + "\\n")
event("start")
time.sleep(0.8)
if os.environ.get("FAIL_ACCESSION") == accession:
    event("end")
    sys.exit(1)
(output / "mash" / "mash.dist").write_text(
    "OTHER:21\\tquery\\t0.22\\t0\\t50/10000\\n"
)
event("end")
"""
        )
        hobrac.chmod(0o755)
        mash = bin_dir / "mash"
        mash.write_text("#!/bin/sh\necho 2.3\n")
        mash.chmod(0o755)
        self.env = {
            **os.environ,
            "PATH": f"{bin_dir}:{Path(sys.executable).parent}:{os.environ['PATH']}",
            "TEST_EVENTS": str(self.events),
        }

    def workflow(self, *extra, expect_success=True):
        result = subprocess.run(
            [
                "bash",
                str(self.scripts / "workflow.sh"),
                str(self.config),
                "--executor",
                "local",
                "--cores",
                "2",
                *extra,
            ],
            env=self.env,
            capture_output=True,
            text=True,
            timeout=120,
        )
        if expect_success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def test_pipeline_cleanup_and_noop_restart(self):
        self.workflow()
        events = [json.loads(line) for line in self.events.read_text().splitlines()]
        starts = {
            (e["kind"], e["accession"]): e["time"]
            for e in events
            if e["event"] == "start"
        }
        ends = {
            (e["kind"], e["accession"]): e["time"]
            for e in events
            if e["event"] == "end"
        }
        for row in self.rows:
            accession = row["accession"]
            self.assertGreaterEqual(
                starts["mash", accession], ends["download", accession]
            )
            run = self.work / "runs" / accession
            self.assertFalse((run / f"{accession}.fna.gz").exists())
            self.assertFalse((run / "hobrac/assembly" / f"{accession}.fna").exists())
            self.assertTrue(
                (run / "hobrac/assembly" / f"{accession}.chr_rename.tsv").exists()
            )
            self.assertTrue((run / f"{accession}.fna.gz.download.json").exists())
            self.assertTrue((run / "completed.json").exists())
            self.assertTrue((run / "execution.json").exists())
            self.assertTrue((run / "hobrac/mash/mash_db.msh").resolve().exists())
        self.assertTrue(
            any(
                starts["download", row["accession"]]
                < starts["mash", other["accession"]]
                < ends["download", row["accession"]]
                for row in self.rows
                for other in self.rows
            ),
            "Mash should overlap another genome's download, not wait for all downloads",
        )
        active = peak = 0
        for event in sorted(events, key=lambda e: e["time"]):
            active += 1 if event["event"] == "start" else -1
            peak = max(peak, active)
        self.assertLessEqual(peak, 2)
        self.assertEqual(len(starts), 8)
        self.assertEqual(
            len(survey.read_tsv(self.work / "report/nearest_references.tsv")), 4
        )
        previous = self.events.read_text()
        self.workflow()
        self.assertEqual(
            self.events.read_text(),
            previous,
            "Restart must not redownload deleted temp inputs",
        )

    def test_failed_genome_is_retained_and_retry_only_reruns_it(self):
        failed = self.rows[0]["accession"]
        self.env["FAIL_ACCESSION"] = failed
        self.workflow(expect_success=False)
        run = self.work / "runs" / failed
        self.assertTrue((run / f"{failed}.fna.gz").exists())
        self.assertTrue((run / "hobrac/assembly" / f"{failed}.fna").exists())
        self.assertFalse((run / "completed.json").exists())
        results = survey.read_tsv(self.work / "report/nearest_references.tsv")
        self.assertEqual([row["status"] for row in results].count("failed"), 1)
        self.assertEqual((self.work / "report/retry_task_ids.txt").read_text(), "0\n")
        before = self.events.read_text().splitlines()
        del self.env["FAIL_ACCESSION"]
        self.workflow()
        added = [
            json.loads(line)
            for line in self.events.read_text().splitlines()[len(before) :]
        ]
        self.assertEqual(
            [(e["kind"], e["accession"]) for e in added if e["event"] == "start"],
            [("mash", failed)],
        )
        self.assertFalse((run / f"{failed}.fna.gz").exists())
        self.assertEqual((self.work / "report/retry_task_ids.txt").read_text(), "")

    def test_legacy_success_is_imported_without_recomputing_mash(self):
        self.workflow()
        # Older surveys had status.json and mash.dist, but no Snakemake markers.
        for row in self.rows:
            (self.work / "runs" / row["accession"] / "completed.json").unlink()
        shutil.rmtree(self.work / ".snakemake")
        before = self.events.read_text().splitlines()
        self.workflow()
        added = [
            json.loads(line)
            for line in self.events.read_text().splitlines()[len(before) :]
        ]
        self.assertFalse(any(event["kind"] == "mash" for event in added))
        for row in self.rows:
            self.assertTrue(
                (self.work / "runs" / row["accession"] / "completed.json").exists()
            )

    def test_thousand_genome_slurm_dag_dry_run(self):
        rows = [
            {
                **self.rows[i % len(self.rows)],
                "task_id": str(i),
                "accession": f"GCA_{i + 100:09}.1",
            }
            for i in range(1000)
        ]
        survey.write_tsv(self.work / "manifest.tsv", rows, ("task_id", *survey.FIELDS))
        snapshot = json.loads((self.work / "snapshot.json").read_text())
        snapshot["count"] = 1000
        snapshot["files"]["manifest.tsv"] = survey.sha256(self.work / "manifest.tsv")
        survey.save_json(self.work / "snapshot.json", snapshot)
        result = self.workflow("--executor", "slurm", "--dry-run")
        output = result.stdout + result.stderr
        for rule, count in (("download_genome", 1000), ("mash", 1000), ("report", 1)):
            self.assertRegex(output, rf"(?m)^{rule}\s+{count}\s*$")
        self.assertFalse(self.events.exists())
        self.assertFalse((self.work / "runs").exists())

    def test_cleanup_can_be_disabled(self):
        self.workflow("--config", "cleanup=false")
        for row in self.rows:
            run = self.work / "runs" / row["accession"]
            self.assertTrue((run / f"{row['accession']}.fna.gz").exists())
            self.assertTrue(
                (run / "hobrac/assembly" / f"{row['accession']}.fna").exists()
            )


if __name__ == "__main__":
    unittest.main()
