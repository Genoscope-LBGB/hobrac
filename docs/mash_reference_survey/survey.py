#!/usr/bin/env python3
"""Frozen, species-balanced NCBI query sampling and a Mash-only hobrac survey."""

import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import random
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import Request, urlopen

GROUPS = (
    "invertebrate",
    "vertebrate_mammalian",
    "vertebrate_other",
    "plant",
    "fungi",
    "protozoa",
)
MASH_URL = "https://www.genoscope.cns.fr/lbgb/mash/"
RANKS = ("species", "genus", "family", "order", "class", "phylum")
FIELDS = (
    "accession",
    "taxid",
    "species_taxid",
    "organism_name",
    "phylum",
    "assembly_level",
    "seq_rel_date",
    "fasta_url",
)


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_json(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def write_tsv(path, rows, fields):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def download(url, path):
    """Atomic downloads with bounded retries and content-addressed receipts."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    receipt = path.with_suffix(path.suffix + ".download.json")
    if path.exists():
        if not receipt.exists():
            raise ValueError(f"Untracked download: {path}; use a fresh work directory")
        recorded = json.loads(receipt.read_text())
        if recorded["url"] != url or recorded["sha256"] != sha256(path):
            raise ValueError(f"Download provenance mismatch: {path}")
        return
    tmp = path.with_suffix(path.suffix + ".part")
    for attempt in range(3):
        try:
            request = Request(url, headers={"User-Agent": "hobrac-mash-survey/1"})
            with urlopen(request, timeout=120) as response, tmp.open("wb") as out:
                shutil.copyfileobj(response, out)
                length = response.headers.get("Content-Length")
                headers = {
                    key: response.headers.get(key) for key in ("ETag", "Last-Modified")
                }
            if tmp.stat().st_size == 0:
                raise ValueError(f"Empty download: {url}")
            if length and tmp.stat().st_size != int(length):
                raise ValueError(f"Truncated download: {url}")
            digest = sha256(tmp)
            save_json(
                receipt,
                {
                    "url": url,
                    "retrieved_utc": now(),
                    "sha256": digest,
                    "bytes": tmp.stat().st_size,
                    "headers": headers,
                },
            )
            tmp.replace(path)
            return
        except (OSError, ValueError):
            tmp.unlink(missing_ok=True)
            if attempt == 2:
                raise
            time.sleep(5 * (attempt + 1))


class Taxonomy:
    """NCBI ranked ancestry, including merged taxids; no live taxonomy queries."""

    def __init__(self, directory):
        directory = Path(directory)
        self.nodes = {}
        with (directory / "nodes.dmp").open() as handle:
            for line in handle:
                fields = [field.strip() for field in line.split("|")]
                self.nodes[int(fields[0])] = (int(fields[1]), sys.intern(fields[2]))
        self.merged = {}
        with (directory / "merged.dmp").open() as handle:
            for line in handle:
                old, new, *_ = line.split("|")
                self.merged[int(old)] = int(new)
        phyla = {taxid for taxid, (_, rank) in self.nodes.items() if rank == "phylum"}
        self.phylum_names = {}
        with (directory / "names.dmp").open() as handle:
            for line in handle:
                fields = [field.strip() for field in line.split("|")]
                taxid = int(fields[0])
                if taxid in phyla and fields[3] == "scientific name":
                    self.phylum_names[taxid] = fields[1]

    def resolve(self, taxid):
        taxid = int(taxid)
        seen = set()
        while taxid in self.merged and taxid not in seen:
            seen.add(taxid)
            taxid = self.merged[taxid]
        return taxid

    @lru_cache(maxsize=200000)
    def lineage(self, taxid):
        taxid = self.resolve(taxid)
        result, ancestors = {}, set()
        while taxid in self.nodes and taxid not in ancestors:
            ancestors.add(taxid)
            parent, rank = self.nodes[taxid]
            result.setdefault(rank, taxid)
            taxid = parent
        return result, frozenset(ancestors)


def assembly_rows(path):
    with Path(path).open() as handle:
        header = None
        for line in handle:
            if line.startswith("# assembly_accession\t"):
                header = line[2:].rstrip("\n").split("\t")
            elif not line.startswith("#") and line.strip():
                if header is None:
                    raise ValueError(f"Missing NCBI assembly-summary header in {path}")
                yield dict(zip(header, line.rstrip("\n").split("\t")))
        if header is None:
            raise ValueError(f"Not an NCBI assembly summary: {path}")


def balanced_sample(rows, count, seed):
    """Water-fill quotas; seed randomises remainder allocation and species draws."""
    pools = defaultdict(list)
    for row in sorted(rows, key=lambda r: (r["phylum"], r["accession"])):
        pools[row["phylum"]].append(row)
    if count <= 0 or count > len(rows):
        raise ValueError(
            f"Requested {count} queries, but only {len(rows)} eligible species"
        )
    rng = random.Random(seed)
    phyla = sorted(pools)
    rng.shuffle(phyla)
    for phylum in phyla:
        rng.shuffle(pools[phylum])
    selected = []
    while len(selected) < count:
        for phylum in phyla:
            if pools[phylum] and len(selected) < count:
                selected.append(pools[phylum].pop())
    return sorted(selected, key=lambda r: (r["phylum"], r["accession"]))


def assembly_priority(row):
    # One deterministic representative per species, before sampling species.
    return (
        {"Complete Genome": 0, "Chromosome": 1}[row["assembly_level"]],
        {"reference genome": 0, "representative genome": 1}.get(
            row.get("refseq_category"), 2
        ),
        -int(row["seq_rel_date"].replace("/", "").replace("-", "")),
        row["assembly_accession"],
    )


def prepare(args):
    work = args.work
    if (work / "snapshot.json").exists():
        raise ValueError(
            "Already prepared; reuse this snapshot or choose a NEW work directory"
        )
    metadata = work / "metadata"
    metadata.mkdir(parents=True, exist_ok=True)
    download(urljoin(MASH_URL, "./"), metadata / "mash_index.html")
    supported = set(
        re.findall(
            r'href="([A-Za-z0-9_-]+)\.msh"',
            (metadata / "mash_index.html").read_text(),
        )
    ) - {"no_returned_phylum"}
    if not supported:
        raise ValueError("No phylum sketches found in the saved Mash index")
    download(
        "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
        metadata / "taxdump.tar.gz",
    )
    taxonomy_dir = metadata / "taxonomy"
    taxonomy_dir.mkdir(exist_ok=True)
    with tarfile.open(metadata / "taxdump.tar.gz", "r:gz") as archive:
        for name in ("nodes.dmp", "names.dmp", "merged.dmp"):
            source = archive.extractfile(name)
            if source is None:
                raise ValueError(f"Missing {name} in taxdump")
            with source, (taxonomy_dir / name).open("wb") as out:
                shutil.copyfileobj(source, out)
    taxonomy = Taxonomy(taxonomy_dir)
    if args.ancestor_taxid not in taxonomy.nodes:
        raise ValueError("Unknown ancestor taxid")
    excluded = Counter()
    representatives = {}
    coverage = defaultdict(Counter)
    for group in GROUPS:
        path = metadata / f"{group}_assembly_summary.txt"
        download(
            f"https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{group}/assembly_summary.txt",
            path,
        )
        for row in assembly_rows(path):
            if row["version_status"] != "latest":
                excluded["not_latest"] += 1
                continue
            if row["assembly_level"] not in ("Chromosome", "Complete Genome"):
                excluded["not_chromosome_or_complete"] += 1
                continue
            if row.get("genome_rep") != "Full" or row.get(
                "excluded_from_refseq", "na"
            ) not in ("", "na"):
                excluded["non_full_or_excluded_from_refseq"] += 1
                continue
            lineage, ancestors = taxonomy.lineage(int(row["taxid"]))
            if args.ancestor_taxid not in ancestors:
                excluded["outside_scope_or_unknown_taxid"] += 1
                continue
            species = lineage.get("species")
            phylum = taxonomy.phylum_names.get(lineage.get("phylum"), "")
            if not species or not phylum:
                excluded["missing_species_or_phylum"] += 1
                continue
            if taxonomy.resolve(row["species_taxid"]) != species:
                excluded["species_metadata_disagreement"] += 1
                continue
            coverage[phylum]["quality_eligible_assemblies"] += 1
            if phylum not in supported:
                excluded["unsupported_phylum"] += 1
                continue
            ftp = row["ftp_path"].replace("ftp://", "https://").rstrip("/")
            if not ftp.startswith("https://ftp.ncbi.nlm.nih.gov/"):
                excluded["missing_download_url"] += 1
                continue
            if not re.fullmatch(r"GCA_\d+\.\d+", row["assembly_accession"]):
                excluded["not_versioned_genbank_accession"] += 1
                continue
            previous = representatives.get(species)
            if previous is None or assembly_priority(row) < assembly_priority(previous):
                representatives[species] = row
    eligible = []
    for species, row in representatives.items():
        lineage, _ = taxonomy.lineage(int(row["taxid"]))
        phylum = taxonomy.phylum_names[lineage["phylum"]]
        ftp = row["ftp_path"].replace("ftp://", "https://").rstrip("/")
        eligible.append(
            {
                "accession": row["assembly_accession"],
                "taxid": row["taxid"],
                "species_taxid": str(species),
                "organism_name": row["organism_name"],
                "phylum": phylum,
                "assembly_level": row["assembly_level"],
                "seq_rel_date": row["seq_rel_date"],
                "fasta_url": f"{ftp}/{ftp.rsplit('/', 1)[1]}_genomic.fna.gz",
            }
        )
        coverage[phylum]["eligible_species"] += 1
    eligible.sort(key=lambda r: (r["phylum"], r["accession"]))
    selected = balanced_sample(eligible, args.count, args.seed)
    for row in selected:
        coverage[row["phylum"]]["sampled_species"] += 1
    for phylum in supported:
        coverage[phylum]  # Include available sketches with zero eligible queries.
    write_tsv(work / "eligible_species.tsv", eligible, FIELDS)
    manifest = [{"task_id": i, **row} for i, row in enumerate(selected)]
    write_tsv(work / "manifest.tsv", manifest, ("task_id", *FIELDS))
    coverage_rows = [
        {
            "phylum": phylum,
            "sketch_available": phylum in supported,
            **{
                key: counts[key]
                for key in (
                    "quality_eligible_assemblies",
                    "eligible_species",
                    "sampled_species",
                )
            },
        }
        for phylum, counts in sorted(coverage.items())
    ]
    write_tsv(
        work / "phylum_coverage.tsv",
        coverage_rows,
        (
            "phylum",
            "sketch_available",
            "quality_eligible_assemblies",
            "eligible_species",
            "sampled_species",
        ),
    )
    save_json(work / "sampling_exclusions.json", dict(excluded))
    for phylum in sorted({row["phylum"] for row in selected}):
        path = work / "sketches" / f"{phylum}.msh"
        download(urljoin(MASH_URL, f"{phylum}.msh"), path)
        info = subprocess.check_output(["mash", "info", str(path)], text=True)
        path.with_suffix(".info.txt").write_text(info)
    tracked = [
        work / "manifest.tsv",
        work / "eligible_species.tsv",
        work / "phylum_coverage.tsv",
        work / "sampling_exclusions.json",
    ]
    tracked += [
        p
        for directory in (metadata, work / "sketches")
        for p in directory.rglob("*")
        if p.is_file()
    ]
    save_json(
        work / "snapshot.json",
        {
            "created_utc": now(),
            "seed": args.seed,
            "count": args.count,
            "ancestor_taxid": args.ancestor_taxid,
            "groups": GROUPS,
            "python": platform.python_version(),
            "script_sha256": sha256(__file__),
            "mash_version": subprocess.check_output(
                ["mash", "--version"], text=True
            ).strip(),
            "files": {str(p.relative_to(work)): sha256(p) for p in sorted(tracked)},
        },
    )
    n_phyla = len({r["phylum"] for r in manifest})
    print(f"Prepared {len(manifest)} queries in {n_phyla} phyla")
    print(f"Review {work / 'phylum_coverage.tsv'} before submitting the array")


def load_snapshot(work, verify_all=False):
    snapshot = json.loads((work / "snapshot.json").read_text())
    paths = snapshot["files"] if verify_all else ["manifest.tsv"]
    for relative in paths:
        if sha256(work / relative) != snapshot["files"][relative]:
            raise ValueError(f"Snapshot checksum mismatch: {relative}")
    return snapshot


def mash_rows(path):
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                raise ValueError(f"Invalid Mash row {number} in {path}")
            reference, query, distance, pvalue, hashes = fields
            accession, separator, taxid = reference.rpartition(":")
            if not separator or not taxid.isdigit():
                raise ValueError(f"Expected accession:taxid in {reference}")
            distance, pvalue = float(distance), float(pvalue)
            shared, total = map(int, hashes.split("/"))
            if not (
                math.isfinite(distance)
                and 0 <= distance <= 1
                and math.isfinite(pvalue)
                and 0 <= pvalue <= 1
                and 0 <= shared <= total
                and total > 0
            ):
                raise ValueError(f"Invalid Mash values in row {number}")
            yield {
                "reference_accession": accession,
                "reference_taxid": taxid,
                "distance": distance,
                "pvalue": pvalue,
                "shared_hashes": shared,
                "sketch_hashes": total,
            }


def download_genome(url, path):
    """Quarantine invalid responses; retry the same accession, never resample."""
    for attempt in range(3):
        download(url, path)
        try:
            with gzip.open(path, "rb") as handle:
                if handle.read(1) != b">":
                    raise ValueError("Downloaded input is not a gzipped FASTA")
                for _ in iter(lambda: handle.read(1024 * 1024), b""):
                    pass
            return
        except (OSError, EOFError, ValueError):
            rejected = path.parent / "rejected_downloads" / str(time.time_ns())
            rejected.mkdir(parents=True)
            path.replace(rejected / path.name)
            receipt = path.with_suffix(path.suffix + ".download.json")
            if receipt.exists():
                receipt.replace(rejected / receipt.name)
            if attempt == 2:
                raise
            time.sleep(5)


def container_environment(container, work):
    env = dict(os.environ)
    if container != "none":
        # Preserve site binds while explicitly exposing the shared sketch target.
        key = f"{container.upper()}_BIND"
        inherited = env.get(key) or env.get(key + "PATH", "")
        env[key] = ",".join(part for part in (inherited, f"{work}:{work}") if part)
    return env


def run_one(args):
    work = args.work
    snapshot = load_snapshot(work)
    manifest = read_tsv(work / "manifest.tsv")
    if not 0 <= args.task_id < len(manifest):
        raise ValueError("Task ID outside manifest range")
    row = manifest[args.task_id]
    run = work / "runs" / row["accession"]
    run.mkdir(parents=True, exist_ok=True)
    # Prevent two array jobs from writing the same query simultaneously.
    import fcntl

    with (run / "task.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        state = run / "status.json"
        if state.exists() and json.loads(state.read_text())["status"] == "success":
            recorded = json.loads(state.read_text())
            result = run / "hobrac" / "mash" / "mash.dist"
            if sha256(result) != recorded["mash_dist_sha256"]:
                raise ValueError(
                    "Completed Mash result changed; use a fresh run directory"
                )
            print(f"Already complete: {row['accession']}")
            return
        save_json(
            state,
            {
                "status": "running",
                "started_utc": now(),
                "task_id": args.task_id,
            },
        )
        try:
            sketch = work / "sketches" / f"{row['phylum']}.msh"
            if sha256(sketch) != snapshot["files"][str(sketch.relative_to(work))]:
                raise ValueError("Frozen reference sketch checksum mismatch")
            genome = run / f"{row['accession']}.fna.gz"
            download_genome(row["fasta_url"], genome)
            output = run / "hobrac"
            mash_dir = output / "mash"
            mash_dir.mkdir(parents=True, exist_ok=True)
            target = mash_dir / "mash_db.msh"
            if not target.exists():
                target.symlink_to(sketch)
            if target.resolve() != sketch:
                raise ValueError("Run points to a different reference sketch")
            (mash_dir / "mash_db.msh.phylum").write_text(row["phylum"] + "\n")
            name = re.sub(r"[^A-Za-z0-9 ._-]", "_", row["organism_name"])
            command = [
                "hobrac",
                "-a",
                str(genome),
                "-n",
                name,
                "-t",
                row["taxid"],
                "-o",
                str(output),
                "-e",
                "local",
                "--stop-after-mash",
                "--rerun-incomplete",
            ]
            if args.container != "none":
                command.append(f"--use-{args.container}")
            versions = {}
            for package in ("HoBRAC", "snakemake"):
                try:
                    versions[package] = importlib.metadata.version(package)
                except importlib.metadata.PackageNotFoundError:
                    versions[package] = "unavailable; record installation commit"
            mash_version = subprocess.check_output(
                ["mash", "--version"], text=True
            ).strip()
            save_json(
                run / "execution.json",
                {
                    "command": command,
                    "command_shell": shlex.join(command),
                    "python": platform.python_version(),
                    "packages": versions,
                    "mash_version": mash_version,
                    "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
                    "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
                    "script_sha256": sha256(__file__),
                    "started_utc": now(),
                },
            )
            with (run / "hobrac.log").open("a") as log:
                log.write(f"\n{now()} {shlex.join(command)}\n")
                log.flush()
                subprocess.run(
                    command,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                    env=container_environment(args.container, work),
                )
            result = mash_dir / "mash.dist"
            if sum(1 for _ in mash_rows(result)) == 0:
                raise ValueError("Mash produced no comparisons")
            # Detect both mutation of the snapshot and replacement of its symlink.
            expected_sketch_hash = snapshot["files"][str(sketch.relative_to(work))]
            if target.resolve() != sketch or sha256(sketch) != expected_sketch_hash:
                raise ValueError("Reference sketch changed during execution")
            save_json(
                state,
                {
                    "status": "success",
                    "finished_utc": now(),
                    "task_id": args.task_id,
                    "mash_dist_sha256": sha256(result),
                    "query_sha256": sha256(genome),
                },
            )
        except Exception as exc:
            save_json(
                state,
                {
                    "status": "failed",
                    "finished_utc": now(),
                    "task_id": args.task_id,
                    "error": str(exc),
                },
            )
            raise


def closest(rows):
    return min(
        rows, key=lambda r: (r["distance"], r["reference_accession"]), default=None
    )


def select_nearest(rows, query, taxonomy):
    """Separate native hobrac eligibility from stricter species-level exclusion."""
    native, strict = [], []
    seen = set()
    counts = Counter()
    species = int(query["species_taxid"])
    for row in rows:
        counts["reference_rows"] += 1
        if row["reference_taxid"] == query["taxid"] or row["distance"] == 0:
            counts["native_excluded"] += 1
            continue
        if row["reference_accession"] in seen:
            counts["duplicate_accessions"] += 1
            continue
        seen.add(row["reference_accession"])
        native.append(row)
        lineage, _ = taxonomy.lineage(int(row["reference_taxid"]))
        ref_species = lineage.get("species")
        if ref_species is None:
            counts["unresolved_reference_species"] += 1
        elif row["reference_accession"] == query["accession"] or ref_species == species:
            counts["additional_same_species_excluded"] += 1
        else:
            strict.append(row)
    counts["eligible_references"] = len(strict)
    # Native ties preserve input order; strict survey ties use accession.
    native_best = min(native, key=lambda r: r["distance"], default=None)
    return native_best, closest(strict), counts


def relationship(query_taxid, reference_taxid, taxonomy):
    query, _ = taxonomy.lineage(int(query_taxid))
    reference, _ = taxonomy.lineage(int(reference_taxid))
    for rank in RANKS:
        if rank in query and query[rank] == reference.get(rank):
            return f"same_{rank}"
    return "unresolved_or_different_phylum"


def summarize(args):
    work = args.work
    load_snapshot(work, verify_all=True)
    taxonomy = Taxonomy(work / "metadata" / "taxonomy")
    manifest = read_tsv(work / "manifest.tsv")
    results = []
    metrics = (
        "reference_accession",
        "reference_taxid",
        "distance",
        "pvalue",
        "shared_hashes",
        "sketch_hashes",
    )
    for query in manifest:
        row = dict(query)
        run = work / "runs" / query["accession"]
        state = run / "status.json"
        row["status"] = "not_run"
        row["error"] = ""
        if state.exists():
            status = json.loads(state.read_text())
            row["status"] = status["status"]
            row["error"] = status.get("error", "")
            if status["status"] == "success":
                try:
                    path = run / "hobrac" / "mash" / "mash.dist"
                    if sha256(path) != status["mash_dist_sha256"]:
                        raise ValueError("Mash output checksum mismatch")
                    native, strict, counts = select_nearest(
                        mash_rows(path), query, taxonomy
                    )
                    row.update(counts)
                    row["hobrac_nearest_accession"] = (
                        native["reference_accession"] if native else ""
                    )
                    row["hobrac_nearest_distance"] = (
                        native["distance"] if native else ""
                    )
                    if strict is None:
                        row["status"] = "no_eligible_reference"
                    else:
                        row.update(strict)
                        row["relationship"] = relationship(
                            query["taxid"], strict["reference_taxid"], taxonomy
                        )
                        row["signal"] = (
                            "no_shared_hashes"
                            if strict["shared_hashes"] == 0
                            else "pvalue_gt_0.01"
                            if strict["pvalue"] > 0.01
                            else "shared_hashes_detected"
                        )
                except (OSError, ValueError) as exc:
                    row["status"], row["error"] = "invalid_result", str(exc)
        results.append(row)
    fields = (
        "task_id",
        *FIELDS,
        "status",
        "error",
        *metrics,
        "relationship",
        "signal",
        "reference_rows",
        "native_excluded",
        "duplicate_accessions",
        "additional_same_species_excluded",
        "unresolved_reference_species",
        "eligible_references",
        "hobrac_nearest_accession",
        "hobrac_nearest_distance",
    )
    report = work / "report"
    report.mkdir(exist_ok=True)
    write_tsv(report / "nearest_references.tsv", results, fields)
    failed = [
        r["task_id"] for r in results if r["status"] in ("not_run", "failed", "running")
    ]
    reset = [r["task_id"] for r in results if r["status"] == "invalid_result"]
    for name, ids in (
        ("retry_task_ids.txt", failed),
        ("reset_required_task_ids.txt", reset),
    ):
        (report / name).write_text(",".join(ids) + ("\n" if ids else ""))
    counts = Counter(r["status"] for r in results)
    save_json(
        report / "summary.json",
        {
            "created_utc": now(),
            "statuses": dict(counts),
            "sampled_queries": len(results),
            "script_sha256": sha256(__file__),
            "snapshot_sha256": sha256(work / "snapshot.json"),
        },
    )
    plot_results(results, report)
    print(json.dumps(dict(counts), indent=2))
    print(f"Tables and plots: {report}")


def plot_results(results, report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    phyla = sorted({r["phylum"] for r in results})
    fig, ax = plt.subplots(figsize=(11, max(4, 0.42 * len(phyla) + 1.5)))
    rng = random.Random(0)
    labels, summaries = [], []
    for position, phylum in enumerate(phyla, 1):
        group = [r for r in results if r["phylum"] == phylum]
        valid = [r for r in group if r["status"] == "success"]
        detected = [float(r["distance"]) for r in valid if r["shared_hashes"] > 0]
        if len(detected) >= 2:
            ax.boxplot(
                detected,
                positions=[position],
                vert=False,
                widths=0.4,
                showfliers=False,
                manage_ticks=False,
            )
        for row in valid:
            no_signal = row["shared_hashes"] == 0
            color = (
                "#c53d3d"
                if no_signal
                else ("#d28a00" if row["pvalue"] > 0.01 else "#266da8")
            )
            ax.scatter(
                row["distance"],
                position + rng.uniform(-0.12, 0.12),
                color=color,
                marker="x" if no_signal else "o",
                s=27,
                zorder=3,
            )
        labels.append(f"{phylum} ({len(valid)}/{len(group)})")
        statuses = Counter(r["status"] for r in group)
        summaries.append(
            {
                "phylum": phylum,
                "sampled": len(group),
                "with_eligible_reference": len(valid),
                "no_shared_hashes": sum(r["shared_hashes"] == 0 for r in valid),
                "pvalue_gt_0.01": sum(r["pvalue"] > 0.01 for r in valid),
                "no_eligible_reference": statuses["no_eligible_reference"],
                "incomplete_or_failed": sum(
                    v
                    for k, v in statuses.items()
                    if k not in ("success", "no_eligible_reference")
                ),
            }
        )
    ax.set_yticks(range(1, len(phyla) + 1), labels)
    ax.invert_yaxis()
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel(
        "Nearest reference Mash distance\n"
        "Native hobrac eligibility + same-species exclusion (zero distances excluded)"
    )
    ax.set_title(
        "Reference availability by query phylum\n"
        "Labels: queries with an eligible reference / sampled"
    )
    ax.grid(axis="x", alpha=0.2)
    ax.legend(
        handles=[
            Line2D(
                [],
                [],
                color="#266da8",
                marker="o",
                linestyle="",
                label="Shared hashes; P ≤ 0.01",
            ),
            Line2D(
                [],
                [],
                color="#d28a00",
                marker="o",
                linestyle="",
                label="Shared hashes; P > 0.01",
            ),
            Line2D(
                [],
                [],
                color="#c53d3d",
                marker="x",
                linestyle="",
                label="No shared hashes (D = 1 sentinel)",
            ),
        ],
        loc="lower right",
        fontsize=8,
    )
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(report / f"distance_by_phylum.{suffix}", dpi=200)
    plt.close(fig)
    write_tsv(report / "phylum_summary.tsv", summaries, tuple(summaries[0]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser(
        "prepare", help="Freeze metadata, sample species, fetch sketches"
    )
    prep.add_argument("--count", type=int, default=100)
    prep.add_argument("--seed", type=int, default=20260717)
    prep.add_argument(
        "--ancestor-taxid",
        type=int,
        default=2759,
        help="Scope: 2759=Eukaryota, 33208=Metazoa",
    )
    prep.set_defaults(func=prepare)
    run = sub.add_parser("run", help="Download and run one zero-based manifest task")
    run.add_argument("--task-id", type=int, required=True)
    run.add_argument(
        "--container", choices=("none", "apptainer", "singularity"), default="none"
    )
    run.set_defaults(func=run_one)
    summary = sub.add_parser(
        "summarize", help="Species-level exclusions, tables and plots"
    )
    summary.set_defaults(func=summarize)
    for command in (prep, run, summary):
        command.add_argument(
            "--work", type=lambda p: Path(p).expanduser().resolve(), required=True
        )
    args = parser.parse_args()
    # hobrac itself interpolates paths into shell commands without consistent quoting.
    if not re.fullmatch(r"[A-Za-z0-9_./-]+", str(args.work)):
        parser.error(
            "Use an absolute work path containing only letters, digits, _, -, . and /"
        )
    try:
        args.func(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
