"""Offline tests: no NCBI downloads, Slurm jobs, or real hobrac execution."""

import argparse
import gzip
import io
import json
import os
import subprocess
import tarfile
import tempfile
import unittest
from collections import Counter
from pathlib import Path
from unittest.mock import patch

import survey

NODES = [
    (1, 1, "no rank"),
    (2759, 1, "superkingdom"),
    (10, 2759, "phylum"),
    (100, 10, "genus"),
    (11, 100, "species"),
    (111, 11, "subspecies"),
    (12, 100, "species"),
    (20, 2759, "phylum"),
    (200, 20, "genus"),
    (21, 200, "species"),
    (22, 200, "species"),
    (30, 2759, "phylum"),
    (31, 30, "species"),
]
NODE_TEXT = "".join(
    f"{taxid}\t|\t{parent}\t|\t{rank}\t|\n" for taxid, parent, rank in NODES
)
NAME_TEXT = "".join(
    f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n"
    for taxid, name in ((10, "Mollusca"), (20, "Arthropoda"), (30, "Unsupported"))
)
MERGED_TEXT = "999\t|\t111\t|\n"
HEADER = (
    "assembly_accession",
    "refseq_category",
    "taxid",
    "species_taxid",
    "organism_name",
    "version_status",
    "assembly_level",
    "genome_rep",
    "seq_rel_date",
    "ftp_path",
    "excluded_from_refseq",
)


def metadata_row(accession, taxid, species=None, level="Chromosome", date="2025/01/01"):
    return dict(
        zip(
            HEADER,
            (
                accession,
                "na",
                str(taxid),
                str(species or taxid),
                f"Test species {taxid}",
                "latest",
                level,
                "Full",
                date,
                f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{accession}_Test",
                "na",
            ),
        )
    )


def reference(accession, taxid, distance, shared=50):
    return {
        "reference_accession": accession,
        "reference_taxid": str(taxid),
        "distance": distance,
        "pvalue": 0.0 if shared else 1.0,
        "shared_hashes": shared,
        "sketch_hashes": 10000,
    }


class Response(io.BytesIO):
    def __init__(self, data, length):
        super().__init__(data)
        self.headers = {"Content-Length": str(length)}


class SurveyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.work = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)

    def taxonomy(self):
        directory = self.work / "taxonomy"
        directory.mkdir(exist_ok=True)
        for name, text in (
            ("nodes.dmp", NODE_TEXT),
            ("names.dmp", NAME_TEXT),
            ("merged.dmp", MERGED_TEXT),
        ):
            (directory / name).write_text(text)
        return survey.Taxonomy(directory)

    def fake_download(self, url, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.name == "taxdump.tar.gz":
            with tarfile.open(path, "w:gz") as archive:
                for name, text in (
                    ("nodes.dmp", NODE_TEXT),
                    ("names.dmp", NAME_TEXT),
                    ("merged.dmp", MERGED_TEXT),
                ):
                    data = text.encode()
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    archive.addfile(info, io.BytesIO(data))
        elif path.name == "mash_index.html":
            path.write_text(
                "".join(
                    f'<a href="{name}.msh">{name}</a>'
                    for name in (
                        "Mollusca",
                        "Arthropoda",
                        "EmptyPhylum",
                        "no_returned_phylum",
                    )
                )
            )
        elif path.name.endswith("assembly_summary.txt"):
            rows = []
            if path.name.startswith("invertebrate_"):
                rows = [
                    metadata_row(f"GCA_{i:09}.1", taxid)
                    for i, taxid in enumerate((11, 12, 21, 22, 31), 1)
                ]
                # Same species, newer version, subspecies taxid: prefer this assembly.
                rows.append(metadata_row("GCA_000000001.2", 111, 11, date="2026/01/01"))
                rows.append(metadata_row("GCA_000000006.1", 12, level="Scaffold"))
                bad = metadata_row("GCA_000000007.1", 12)
                bad["excluded_from_refseq"] = "contaminated"
                rows.append(bad)
            path.write_text(
                "#assembly_accession\t"
                + "\t".join(HEADER[1:])
                + "\n"
                + "".join("\t".join(row[k] for k in HEADER) + "\n" for row in rows)
            )
        elif path.suffix == ".msh":
            path.write_bytes(b"mock sketch: " + path.name.encode())
        elif path.name.endswith(".fna.gz"):
            with gzip.open(path, "wb") as handle:
                handle.write(b">chr1\nACGTACGTACGT\n")
        else:
            raise AssertionError(f"Unexpected download: {url}")

    def prepare(self, count=4):
        args = argparse.Namespace(
            work=self.work, count=count, seed=17, ancestor_taxid=2759
        )
        with (
            patch.object(survey, "download", side_effect=self.fake_download),
            patch.object(survey.subprocess, "check_output", return_value="2.3"),
        ):
            survey.prepare(args)
        return survey.read_tsv(self.work / "manifest.tsv")

    def test_assembly_rows_accepts_ncbi_header_spacing(self):
        path = self.work / "assembly_summary.txt"
        expected = metadata_row("GCA_000000001.1", 11)
        expected["excluded_from_refseq"] = ""
        for spacing in ("", " ", "  ", "\t"):
            with self.subTest(spacing=spacing):
                path.write_text(
                    "## See README_assembly_summary.txt for column descriptions.\n"
                    + "#"
                    + spacing
                    + "\t".join(HEADER)
                    + "\n"
                    + "\t".join(expected[k] for k in HEADER)
                    + "\n"
                )
                self.assertEqual(list(survey.assembly_rows(path)), [expected])

    def test_assembly_rows_rejects_missing_header(self):
        path = self.work / "assembly_summary.txt"
        for text in (
            "",
            "## Description only\n",
            "#assembly_accession_extra\ttaxid\nGCA_000000001.1\t11\n",
            "GCA_000000001.1\t11\n",
        ):
            with self.subTest(text=text):
                path.write_text(text)
                with self.assertRaisesRegex(ValueError, "assembly.summary"):
                    list(survey.assembly_rows(path))

    def test_balanced_sample_reproducible_and_order_independent(self):
        rows = [
            {"phylum": p, "accession": f"{p}{i}"}
            for p, size in (("A", 10), ("B", 10), ("C", 1))
            for i in range(size)
        ]
        chosen = survey.balanced_sample(rows, 9, 42)
        self.assertEqual(chosen, survey.balanced_sample(list(reversed(rows)), 9, 42))
        self.assertEqual(Counter(r["phylum"] for r in chosen), {"A": 4, "B": 4, "C": 1})
        self.assertNotEqual(chosen, survey.balanced_sample(rows, 9, 100))
        with self.assertRaises(ValueError):
            survey.balanced_sample(rows, 22, 42)
        with self.assertRaises(ValueError):
            survey.balanced_sample(rows, 0, 42)

    def test_taxonomy_subspecies_merged_and_unknown(self):
        taxonomy = self.taxonomy()
        self.assertEqual(taxonomy.lineage(999)[0]["species"], 11)
        self.assertEqual(taxonomy.lineage(111)[0]["phylum"], 10)
        self.assertEqual(taxonomy.lineage(999999)[0], {})
        self.assertEqual(survey.relationship(111, 12, taxonomy), "same_genus")

    def test_species_exclusion_and_unknown_reference_taxonomy(self):
        query = {"accession": "QUERY", "taxid": "11", "species_taxid": "11"}
        rows = [
            reference("QUERY", 11, 0),
            reference("SUBSPECIES", 111, 0.01),
            reference("MERGED", 999, 0.02),
            reference("UNKNOWN", 99999, 0.03),
            reference("ZERO", 12, 0),
            reference("B", 12, 0.2),
            reference("A", 12, 0.2),
            reference("A", 12, 0.3),
        ]
        native, strict, counts = survey.select_nearest(rows, query, self.taxonomy())
        assert native is not None and strict is not None
        self.assertEqual(native["reference_accession"], "SUBSPECIES")
        self.assertEqual(strict["reference_accession"], "A")
        self.assertEqual(counts["additional_same_species_excluded"], 2)
        self.assertEqual(counts["unresolved_reference_species"], 1)
        self.assertEqual(counts["duplicate_accessions"], 1)
        self.assertEqual(counts["native_excluded"], 2)

    def test_no_eligible_reference_is_not_distance_one(self):
        query = {"accession": "QUERY", "taxid": "11", "species_taxid": "11"}
        taxonomy = self.taxonomy()
        _, strict, _ = survey.select_nearest(
            [reference("SELF", 111, 0.01)], query, taxonomy
        )
        self.assertIsNone(strict)
        _, strict, _ = survey.select_nearest(
            [reference("OTHER", 12, 1, 0)], query, taxonomy
        )
        assert strict is not None
        self.assertEqual(strict["distance"], 1)
        self.assertEqual(strict["shared_hashes"], 0)

    def test_mash_parser_rejects_bad_rows(self):
        path = self.work / "mash.dist"
        for text in (
            "broken",
            "REF:12\tquery\tnan\t0\t1/1000\n",
            "REF:12\tquery\t0.2\t0\t1001/1000\n",
            "REF\tquery\t0.2\t0\t1/1000\n",
        ):
            path.write_text(text)
            with self.assertRaises(ValueError):
                list(survey.mash_rows(path))
        path.write_text("REF:12\tquery\t1\t1\t0/10000\n")
        self.assertEqual(list(survey.mash_rows(path))[0]["shared_hashes"], 0)

    def test_prepare_freezes_eligibility_manifest_and_sketches(self):
        rows = self.prepare()
        self.assertEqual(len(rows), 4)
        self.assertEqual(len({r["species_taxid"] for r in rows}), 4)
        self.assertIn("GCA_000000001.2", [r["accession"] for r in rows])
        self.assertNotIn("GCA_000000001.1", [r["accession"] for r in rows])
        coverage = {
            r["phylum"]: r for r in survey.read_tsv(self.work / "phylum_coverage.tsv")
        }
        self.assertEqual(coverage["Unsupported"]["sketch_available"], "False")
        self.assertEqual(coverage["EmptyPhylum"]["sampled_species"], "0")
        survey.load_snapshot(self.work, verify_all=True)
        with self.assertRaises(ValueError):
            self.prepare()
        with (self.work / "manifest.tsv").open("a") as handle:
            handle.write("changed\n")
        with self.assertRaises(ValueError):
            survey.load_snapshot(self.work)

    def test_run_and_summary_keep_missing_queries_visible(self):
        rows = self.prepare()
        query = rows[0]
        other_taxid = 22 if query["species_taxid"] == "21" else 21

        def hobrac(command, **kwargs):
            self.assertEqual(command[command.index("-e") + 1], "local")
            self.assertIn("--stop-after-mash", command)
            output = Path(command[command.index("-o") + 1])
            self.assertTrue((output / "mash" / "mash_db.msh").is_symlink())
            (output / "mash" / "mash.dist").write_text(
                f"SELF:{query['taxid']}\tquery\t0\t0\t10000/10000\n"
                f"OTHER:{other_taxid}\tquery\t0.22\t0\t50/10000\n"
            )

        args = argparse.Namespace(work=self.work, task_id=0, container="none")
        with (
            patch.object(survey, "download", side_effect=self.fake_download),
            patch.object(survey.subprocess, "check_output", return_value="2.3"),
            patch.object(survey.subprocess, "run", side_effect=hobrac) as run,
        ):
            survey.run_one(args)
            survey.run_one(args)
            self.assertEqual(run.call_count, 1)
        with patch.object(survey, "plot_results"):
            survey.summarize(argparse.Namespace(work=self.work))
        results = survey.read_tsv(self.work / "report" / "nearest_references.tsv")
        self.assertEqual(len(results), 4)
        self.assertEqual(results[0]["distance"], "0.22")
        self.assertEqual([r["status"] for r in results], ["success"] + ["not_run"] * 3)
        retry = self.work / "report" / "retry_task_ids.txt"
        self.assertEqual(retry.read_text(), "1,2,3\n")
        # A corrupt successful output must not become a distance observation.
        result = self.work / "runs" / query["accession"] / "hobrac/mash/mash.dist"
        result.write_text("corrupt\n")
        with patch.object(survey, "plot_results"):
            survey.summarize(argparse.Namespace(work=self.work))
        results = survey.read_tsv(self.work / "report" / "nearest_references.tsv")
        self.assertEqual(results[0]["status"], "invalid_result")
        self.assertEqual(retry.read_text(), "1,2,3\n")
        reset = self.work / "report" / "reset_required_task_ids.txt"
        self.assertEqual(reset.read_text(), "0\n")

    def test_container_binds_preserve_site_configuration(self):
        for container in ("apptainer", "singularity"):
            key = container.upper() + "_BIND"
            with patch.dict(os.environ, {key: "/site:/site"}):
                env = survey.container_environment(container, self.work)
                self.assertEqual(env[key], f"/site:/site,{self.work}:{self.work}")
        with patch.dict(os.environ, {"APPTAINER_BINDPATH": "/site:/site"}, clear=True):
            env = survey.container_environment("apptainer", self.work)
            self.assertEqual(
                env["APPTAINER_BIND"], f"/site:/site,{self.work}:{self.work}"
            )

    def test_invalid_genome_is_quarantined_and_redownloaded(self):
        path = self.work / "query.fna.gz"
        attempts = []

        def fetch(url, target):
            attempts.append(url)
            if len(attempts) == 1:
                target.write_bytes(b"not a gzipped FASTA")
            else:
                self.fake_download(url, target)

        with (
            patch.object(survey, "download", side_effect=fetch),
            patch.object(survey.time, "sleep"),
        ):
            survey.download_genome("https://example.org/query.fna.gz", path)
        self.assertEqual(len(attempts), 2)
        self.assertEqual(attempts[0], attempts[1])
        self.assertEqual(len(list((self.work / "rejected_downloads").glob("*/*"))), 1)
        with gzip.open(path, "rb") as handle:
            self.assertTrue(handle.read().startswith(b">chr1"))

    def test_run_records_failure(self):
        self.prepare()
        args = argparse.Namespace(work=self.work, task_id=0, container="none")
        with patch.object(survey, "download", side_effect=OSError("offline")):
            with self.assertRaises(OSError):
                survey.run_one(args)
        query = survey.read_tsv(self.work / "manifest.tsv")[0]
        state_path = self.work / "runs" / query["accession"] / "status.json"
        state = json.loads(state_path.read_text())
        self.assertEqual(state["status"], "failed")
        self.assertEqual(state["error"], "offline")

    def test_download_receipt_and_reuse(self):
        response = Response(b"hello", 5)
        path = self.work / "download.txt"
        with patch.object(survey, "urlopen", return_value=response) as fetch:
            survey.download("https://example.org/file", path)
            survey.download("https://example.org/file", path)
            self.assertEqual(fetch.call_count, 1)
        path.write_text("changed")
        with self.assertRaises(ValueError):
            survey.download("https://example.org/file", path)

    def test_truncated_download_never_installed(self):
        def response(*args, **kwargs):
            return Response(b"short", 100)

        path = self.work / "download.txt"
        with (
            patch.object(survey, "urlopen", side_effect=response) as fetch,
            patch.object(survey.time, "sleep"),
        ):
            with self.assertRaises(ValueError):
                survey.download("https://example.org/file", path)
            self.assertEqual(fetch.call_count, 3)
        self.assertFalse(path.exists())
        self.assertFalse(path.with_suffix(".txt.part").exists())

    def test_slurm_array_and_afterany_summary_submission(self):
        self.prepare()
        scripts = Path(survey.__file__).parent
        config = (
            (scripts / "config.example.sh")
            .read_text()
            .replace("/path/to/shared/scratch/hobrac_mash_survey", str(self.work))
        )
        config_path = self.work / "config.sh"
        config_path.write_text(config)
        bin_dir = self.work / "bin"
        bin_dir.mkdir()
        sbatch = bin_dir / "sbatch"
        sbatch.write_text(
            '#!/bin/bash\nprintf "%s\\n" "$*" >> "$TEST_LOG"\necho "123;cluster"\n'
        )
        sbatch.chmod(0o755)
        log = self.work / "submissions.log"
        env = {
            **os.environ,
            "PATH": f"{bin_dir}:{os.environ['PATH']}",
            "TEST_LOG": str(log),
        }
        subprocess.run(
            ["bash", str(scripts / "submit.sh"), "array", str(config_path)],
            env=env,
            check=True,
            capture_output=True,
            text=True,
        )
        calls = log.read_text().splitlines()
        self.assertEqual(len(calls), 2)
        self.assertIn("--array=0-3%4", calls[0])
        self.assertIn("--dependency=afterany:123", calls[1])
        self.assertIn(str(scripts / "job.sh"), calls[0])

    def test_plot_includes_empty_phylum_and_no_shared_hashes(self):
        try:
            import matplotlib  # noqa: F401
        except ImportError:
            self.skipTest("Install matplotlib to exercise rendering")
        results = [
            {
                "phylum": "A",
                "status": "success",
                "distance": 0.2,
                "shared_hashes": 30,
                "pvalue": 0.0,
            },
            {
                "phylum": "A",
                "status": "success",
                "distance": 0.3,
                "shared_hashes": 10,
                "pvalue": 0.02,
            },
            {
                "phylum": "A",
                "status": "success",
                "distance": 1.0,
                "shared_hashes": 0,
                "pvalue": 1.0,
            },
            {"phylum": "B", "status": "failed"},
        ]
        survey.plot_results(results, self.work)
        self.assertGreater((self.work / "distance_by_phylum.png").stat().st_size, 1000)
        self.assertGreater((self.work / "distance_by_phylum.pdf").stat().st_size, 1000)
        rows = survey.read_tsv(self.work / "phylum_summary.tsv")
        self.assertEqual(rows[0]["no_shared_hashes"], "1")
        self.assertEqual(rows[1]["incomplete_or_failed"], "1")


if __name__ == "__main__":
    unittest.main()
