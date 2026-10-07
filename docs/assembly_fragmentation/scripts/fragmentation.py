#!/usr/bin/env python3
"""ATLASea assembly-fragmentation experiment: frozen sampling, per-level FASTAs,
full hobrac runs and the runtime/dotplot report."""

import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import json
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
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

BIOPROJECT = "PRJEB64126"  # ATLASea umbrella BioProject
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/"
FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
TAXDUMP = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
USER_AGENT = "hobrac-fragmentation/1"
LEVELS = ("contigs", "scaffolds", "chromosomes")
# N runs at least this long are assembly gaps; shorter runs stay inside contigs.
MIN_GAP = 10
ORGANELLES = {
    "Apicoplast",
    "Chloroplast",
    "Chromatophore",
    "Cyanelle",
    "Kinetoplast",
    "Mitochondrion",
    "Nucleomorph",
    "Plastid",
}
SEQUENCE_ROLES = {
    "assembled-molecule": "chromosome",
    "unlocalized-scaffold": "unlocalized",
    "unplaced-scaffold": "unplaced",
}
# Options that would change what is being measured: the input, the automatic
# reference choice, the frozen MASH release, BUSCO reruns or the minimap2 step.
FORBIDDEN_HOBRAC_ARGS = {
    "-a",
    "--assembly",
    "-n",
    "--name",
    "-t",
    "--taxid",
    "-o",
    "--output",
    "-r",
    "--reference",
    "-e",
    "--executor",
    "--profile",
    "--mash-db-version",
    "--busco-assembly",
    "--busco-reference",
    "--stop-after-mash",
    "--skip-genomic",
    "--ref-count",
}
CANDIDATE_FIELDS = (
    "accession",
    "assembly_name",
    "organism_name",
    "taxid",
    "species_taxid",
    "phylum",
    "assembly_level",
    "assembly_type",
    "release_date",
    "refseq_category",
    "total_length",
    "contig_n50",
    "scaffold_n50",
    "number_of_contigs",
    "number_of_scaffolds",
    "number_of_chromosomes",
    "status",
)
MANIFEST_FIELDS = ("task_id", *CANDIDATE_FIELDS[:-1], "fasta_url")
STEP_GROUPS = (
    "busco_assembly",
    "busco_reference",
    "minimap2",
    "mash",
    "download",
    "other",
)
RUN_FIELDS = (
    "task_id",
    "accession",
    "organism_name",
    "phylum",
    "taxid",
    "level",
    "status",
    "error",
    "attempt",
    "sequences",
    "total_length",
    "n50",
    "n_bases",
    "ncbi_number_of_contigs",
    "hobrac_phylum",
    "reference_accession",
    "reference_distance",
    "reference_pvalue",
    "busco_lineage",
    "busco_complete",
    "busco_single",
    "busco_duplicated",
    "busco_fragmented",
    "busco_missing",
    "busco_markers",
    "busco_anchors",
    "anchor_sequences",
    "significant_pairs",
    "significant_anchor_fraction",
    "hobrac_wall_seconds",
    *(f"{group}_seconds" for group in STEP_GROUPS),
    "node",
    "cpu_model",
    "cpus",
    "slurm_job_id",
    "dotplot_busco",
    "dotplot_busco_bw",
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
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_jsonl(path):
    with Path(path).open() as handle:
        return [json.loads(line) for line in handle if line.strip()]


def http_get(url, method="GET"):
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    if os.environ.get("NCBI_API_KEY"):
        headers["api-key"] = os.environ["NCBI_API_KEY"]
    for attempt in range(3):
        try:
            request = Request(url, headers=headers, method=method)
            with urlopen(request, timeout=120) as response:
                return response.read()
        except OSError:
            if attempt == 2:
                raise
            time.sleep(5 * (attempt + 1))


def download(url, path):
    """Atomic download with bounded retries and a content-addressed receipt."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    receipt = path.with_suffix(path.suffix + ".download.json")
    if path.exists() and receipt.exists():
        recorded = json.loads(receipt.read_text())
        if recorded["url"] == url and recorded["sha256"] == sha256(path):
            return
    tmp = path.with_suffix(path.suffix + ".part")
    for attempt in range(3):
        try:
            request = Request(url, headers={"User-Agent": USER_AGENT})
            with urlopen(request, timeout=120) as response, tmp.open("wb") as out:
                shutil.copyfileobj(response, out)
                length = response.headers.get("Content-Length")
            if tmp.stat().st_size == 0:
                raise ValueError(f"Empty download: {url}")
            if length and tmp.stat().st_size != int(length):
                raise ValueError(f"Truncated download: {url}")
            save_json(
                receipt,
                {
                    "url": url,
                    "retrieved_utc": now(),
                    "sha256": sha256(tmp),
                    "bytes": tmp.stat().st_size,
                },
            )
            tmp.replace(path)
            return
        except (OSError, ValueError):
            tmp.unlink(missing_ok=True)
            if attempt == 2:
                raise
            time.sleep(5 * (attempt + 1))


def download_genome(url, path):
    """Reject anything that is not a complete gzipped FASTA; never resample."""
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
            Path(path).unlink(missing_ok=True)
            if attempt == 2:
                raise
            time.sleep(5)


def api_reports(endpoint):
    """All pages of an NCBI Datasets v2 report, checked against total_count."""
    reports, token = [], None
    while True:
        query = {"page_size": 1000}
        if token:
            query["page_token"] = token
        payload = json.loads(http_get(f"{API}{endpoint}?{urlencode(query)}"))
        reports.extend(payload.get("reports", []))
        token = payload.get("next_page_token")
        if not token:
            break
    expected = payload.get("total_count")
    if expected is not None and int(expected) != len(reports):
        raise ValueError(f"{endpoint}: expected {expected} reports, got {len(reports)}")
    return reports


def fetch_jsonl(endpoint, path):
    reports = api_reports(endpoint)
    if not reports:
        raise ValueError(f"No reports returned by {endpoint}")
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".part")
    with tmp.open("w") as handle:
        for report in reports:
            handle.write(json.dumps(report, sort_keys=True) + "\n")
    tmp.replace(path)
    save_json(
        path.with_suffix(path.suffix + ".download.json"),
        {
            "url": f"{API}{endpoint}",
            "retrieved_utc": now(),
            "sha256": sha256(path),
            "reports": len(reports),
        },
    )
    return reports


class Taxonomy:
    """NCBI ranked ancestry from a frozen taxdump, following merged taxids."""

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
        taxid, seen = int(taxid), set()
        while taxid in self.merged and taxid not in seen:
            seen.add(taxid)
            taxid = self.merged[taxid]
        return taxid

    def lineage(self, taxid):
        taxid = self.resolve(taxid)
        result, seen = {}, set()
        while taxid in self.nodes and taxid not in seen:
            seen.add(taxid)
            parent, rank = self.nodes[taxid]
            result.setdefault(rank, taxid)
            taxid = parent
        return result


def fasta_url(accession, assembly_name):
    digits = accession[4:].split(".")[0]
    name = re.sub(r"[^A-Za-z0-9._-]", "_", assembly_name)
    stem = f"{accession}_{name}"
    return (
        f"{FTP}{accession[:3]}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/"
        f"{stem}/{stem}_genomic.fna.gz"
    )


def candidate_row(report):
    info = report.get("assembly_info", {})
    stats = report.get("assembly_stats", {})
    return {
        "accession": report["accession"],
        "assembly_name": info.get("assembly_name", ""),
        "organism_name": report["organism"]["organism_name"],
        "taxid": str(report["organism"]["tax_id"]),
        "species_taxid": "",
        "phylum": "",
        "assembly_level": info.get("assembly_level", ""),
        "assembly_type": info.get("assembly_type", ""),
        "release_date": info.get("release_date", ""),
        "refseq_category": info.get("refseq_category", ""),
        "total_length": stats.get("total_sequence_length", ""),
        "contig_n50": stats.get("contig_n50", ""),
        "scaffold_n50": stats.get("scaffold_n50", ""),
        "number_of_contigs": stats.get("number_of_contigs", ""),
        "number_of_scaffolds": stats.get("number_of_scaffolds", ""),
        "number_of_chromosomes": stats.get("total_number_of_chromosomes", ""),
        "status": "",
        "_assembly_status": info.get("assembly_status", ""),
    }


def representative_priority(row):
    # One deterministic assembly per species before drawing per phylum.
    return (
        {"Complete Genome": 0, "Chromosome": 1}[row["assembly_level"]],
        0 if row["refseq_category"] == "reference genome" else 1,
        -int(row["release_date"].replace("-", "") or 0),
        row["accession"],
    )


def select_candidates(reports, taxonomy, sketches):
    """Annotate every assembly with its phylum and eligibility status."""
    rows, representatives = [], {}
    for report in reports:
        row = candidate_row(report)
        rows.append(row)
        if not re.fullmatch(r"GCA_\d+\.\d+", row["accession"]):
            row["status"] = "excluded:not_genbank_accession"
            continue
        if row.pop("_assembly_status") != "current":
            row["status"] = "excluded:not_current"
            continue
        if row["assembly_type"] == "alternate-pseudohaplotype":
            row["status"] = "excluded:alternate_haplotype"
            continue
        if row["assembly_level"] not in ("Chromosome", "Complete Genome"):
            row["status"] = "excluded:not_chromosome_level"
            continue
        if taxonomy.resolve(row["taxid"]) not in taxonomy.nodes:
            row["status"] = "excluded:taxid_not_in_taxdump"
            continue
        lineage = taxonomy.lineage(row["taxid"])
        row["phylum"] = taxonomy.phylum_names.get(lineage.get("phylum"), "")
        if not row["phylum"]:
            row["status"] = "excluded:no_phylum"
            continue
        if "species" not in lineage:
            row["status"] = "excluded:no_species_rank"
            continue
        row["species_taxid"] = str(lineage["species"])
        if row["phylum"] not in sketches:
            row["status"] = "excluded:no_mash_sketch"
            continue
        previous = representatives.get(row["species_taxid"])
        if previous is None or representative_priority(row) < representative_priority(
            previous
        ):
            if previous is not None:
                previous["status"] = "excluded:same_species_duplicate"
            representatives[row["species_taxid"]] = row
            row["status"] = "eligible"
        else:
            row["status"] = "excluded:same_species_duplicate"
    for row in rows:
        row.pop("_assembly_status", None)
    return rows


def sample_one_per_phylum(rows, seed):
    """One eligible species per phylum; independent of the API's row order."""
    pools = defaultdict(list)
    for row in sorted(rows, key=lambda r: (r["phylum"], r["accession"])):
        if row["status"] == "eligible":
            pools[row["phylum"]].append(row)
    rng = random.Random(seed)
    return [rng.choice(pools[phylum]) for phylum in sorted(pools)]


def prepare(args):
    from hobrac.mash_database import mash_db_url, resolve_mash_db_version

    work = args.work
    if (work / "snapshot.json").exists():
        raise ValueError(
            "Already prepared; reuse this snapshot or choose a NEW work directory"
        )
    metadata = work / "metadata"
    metadata.mkdir(parents=True, exist_ok=True)

    version = resolve_mash_db_version(args.mash_db_version)
    download(mash_db_url(version), metadata / "mash_index.html")
    sketches = set(
        re.findall(
            r'href="([A-Za-z0-9_-]+)\.msh"', (metadata / "mash_index.html").read_text()
        )
    ) - {"no_returned_phylum"}
    if not sketches:
        raise ValueError(f"No phylum sketches listed for MASH release {version}")

    download(TAXDUMP, metadata / "taxdump.tar.gz")
    taxonomy_dir = metadata / "taxonomy"
    taxonomy_dir.mkdir(exist_ok=True)
    with tarfile.open(metadata / "taxdump.tar.gz", "r:gz") as archive:
        # taxonkit, used by hobrac, needs delnodes.dmp as well.
        for name in ("nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"):
            source = archive.extractfile(name)
            if source is None:
                raise ValueError(f"Missing {name} in taxdump")
            with source, (taxonomy_dir / name).open("wb") as out:
                shutil.copyfileobj(source, out)
    taxonomy = Taxonomy(taxonomy_dir)

    reports = fetch_jsonl(
        f"genome/bioproject/{args.bioproject}/dataset_report",
        metadata / "bioproject_assemblies.jsonl",
    )
    rows = select_candidates(reports, taxonomy, sketches)
    sampled = sample_one_per_phylum(rows, args.seed)
    for row in sampled:
        row["status"] = "sampled"
    rows.sort(key=lambda r: (r["phylum"] or "~", r["accession"]))
    write_tsv(work / "atlasea_assemblies.tsv", rows, CANDIDATE_FIELDS)

    manifest = []
    for task_id, row in enumerate(sorted(sampled, key=lambda r: r["phylum"])):
        url = fasta_url(row["accession"], row["assembly_name"])
        http_get(url, method="HEAD")  # Fail now, not in a Slurm job, on a bad path.
        manifest.append({"task_id": task_id, **row, "fasta_url": url})
    write_tsv(work / "manifest.tsv", manifest, MANIFEST_FIELDS)

    coverage = defaultdict(Counter)
    for row in rows:
        if row["phylum"]:
            coverage[row["phylum"]]["assemblies"] += 1
            if row["status"] in ("eligible", "sampled"):
                coverage[row["phylum"]]["eligible_species"] += 1
    sampled_by_phylum = {row["phylum"]: row["accession"] for row in sampled}
    write_tsv(
        work / "phylum_coverage.tsv",
        [
            {
                "phylum": phylum,
                "assemblies": counts["assemblies"],
                "eligible_species": counts["eligible_species"],
                "mash_sketch": phylum in sketches,
                "sampled_accession": sampled_by_phylum.get(phylum, ""),
            }
            for phylum, counts in sorted(coverage.items())
        ],
        (
            "phylum",
            "assemblies",
            "eligible_species",
            "mash_sketch",
            "sampled_accession",
        ),
    )
    save_json(
        work / "sampling_exclusions.json",
        dict(Counter(r["status"] for r in rows if r["status"].startswith("excluded"))),
    )

    tracked = [
        work / name
        for name in (
            "atlasea_assemblies.tsv",
            "manifest.tsv",
            "phylum_coverage.tsv",
            "sampling_exclusions.json",
        )
    ]
    tracked += sorted(p for p in metadata.rglob("*") if p.is_file())
    try:
        hobrac_version = importlib.metadata.version("HoBRAC")
    except importlib.metadata.PackageNotFoundError:
        hobrac_version = "unavailable; record installation commit"
    save_json(
        work / "snapshot.json",
        {
            "created_utc": now(),
            "bioproject": args.bioproject,
            "seed": args.seed,
            "mash_db_version": version,
            "mash_db_version_requested": args.mash_db_version,
            "min_gap": MIN_GAP,
            "hobrac_version": hobrac_version,
            "python": platform.python_version(),
            "script_sha256": sha256(__file__),
            "files": {str(p.relative_to(work)): sha256(p) for p in tracked},
        },
    )
    print(
        f"Prepared {len(manifest)} genomes ({len(manifest) * len(LEVELS)} hobrac runs) "
        f"with MASH release {version}"
    )
    print(f"Review {work / 'manifest.tsv'} and {work / 'phylum_coverage.tsv'}")


def load_snapshot(work, verify_all=False):
    snapshot = json.loads((work / "snapshot.json").read_text())
    paths = snapshot["files"] if verify_all else ["manifest.tsv"]
    for relative in paths:
        if sha256(work / relative) != snapshot["files"][relative]:
            raise ValueError(f"Snapshot checksum mismatch: {relative}")
    return snapshot


def check_work(args):
    load_snapshot(args.work, verify_all=True)


def task_context(work, task_id):
    snapshot = load_snapshot(work)
    manifest = read_tsv(work / "manifest.tsv")
    if not 0 <= task_id < len(manifest):
        raise ValueError("Task ID outside manifest range")
    row = manifest[task_id]
    run = work / "runs" / row["accession"]
    run.mkdir(parents=True, exist_ok=True)
    return snapshot, row, run


def download_one(args):
    _, row, run = task_context(args.work, args.task_id)
    download_genome(row["fasta_url"], run / f"{row['accession']}.fna.gz")
    fetch_jsonl(
        f"genome/accession/{row['accession']}/sequence_reports",
        run / "sequence_report.jsonl",
    )


def read_fasta(handle):
    name, chunks, seen = None, [], set()
    for line in handle:
        if line.startswith(b">"):
            if name is not None:
                yield name, b"".join(chunks)
            fields = line[1:].split(None, 1)
            if not fields:
                raise ValueError("FASTA record without an identifier")
            name, chunks = fields[0].decode(), []
            if name in seen:
                raise ValueError(f"Duplicate FASTA identifier: {name}")
            seen.add(name)
        elif name is None:
            if line.strip():
                raise ValueError("FASTA does not start with a header")
        else:
            chunks.append(line.rstrip())
    if name is not None:
        yield name, b"".join(chunks)


def write_fasta(handle, name, sequence, width=80):
    lines = b"\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))
    handle.write(b">" + name.encode() + b"\n" + lines + b"\n")


def split_contigs(sequence, min_gap=MIN_GAP):
    return [piece for piece in re.split(rb"[Nn]{%d,}" % min_gap, sequence) if piece]


def n50(lengths):
    total, running = sum(lengths), 0
    for length in sorted(lengths, reverse=True):
        running += length
        if 2 * running >= total:
            return length
    return 0


def classify_sequences(report_rows):
    """FASTA id -> (kind, output name); organelles map to kind 'organelle'."""
    result, used = {}, set()
    for row in report_rows:
        accession = row.get("genbank_accession")
        if not accession:
            raise ValueError(f"Sequence report row without a GenBank accession: {row}")
        if (
            row.get("assembly_unit") == "non-nuclear"
            or row.get("assigned_molecule_location_type") in ORGANELLES
        ):
            result[accession] = ("organelle", accession)
            continue
        kind = SEQUENCE_ROLES.get(row.get("role"))
        if kind is None:
            raise ValueError(
                f"Unsupported sequence role {row.get('role')!r} for {accession}"
            )
        name = accession
        chr_name = row.get("chr_name", "")
        if kind == "chromosome" and chr_name and chr_name.lower() != "na":
            candidate = "chr" + re.sub(r"[^A-Za-z0-9.-]+", "_", chr_name)
            if candidate not in used:
                name = candidate
        if name in used:
            raise ValueError(f"Duplicate sequence name {name}")
        used.add(name)
        result[accession] = (kind, name)
    return result


def build_levels(genome, report_rows, out_dir, min_gap=MIN_GAP):
    """Write the three level FASTAs in one pass and return their statistics."""
    from xopen import xopen

    classes = classify_sequences(report_rows)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp = {level: out_dir / f".{level}.fna.gz.tmp" for level in LEVELS}
    lengths = {level: [] for level in LEVELS}
    n_bases = Counter()
    organelles, seen = 0, set()
    handles = {
        level: xopen(path, "wb", compresslevel=1, format="gz")
        for level, path in tmp.items()
    }
    try:
        with xopen(genome, "rb") as source:
            for accession, sequence in read_fasta(source):
                if accession not in classes:
                    raise ValueError(f"{accession} is absent from the sequence report")
                seen.add(accession)
                kind, name = classes[accession]
                if kind == "organelle":
                    organelles += 1
                    continue
                ns = sequence.count(b"N") + sequence.count(b"n")
                targets = (
                    ("scaffolds", "chromosomes")
                    if kind == "chromosome"
                    else ("scaffolds",)
                )
                for level in targets:
                    write_fasta(handles[level], name, sequence)
                    lengths[level].append(len(sequence))
                    n_bases[level] += ns
                pieces = split_contigs(sequence, min_gap)
                for index, piece in enumerate(pieces, 1):
                    piece_name = name if len(pieces) == 1 else f"{name}_ctg{index}"
                    write_fasta(handles["contigs"], piece_name, piece)
                    lengths["contigs"].append(len(piece))
                    n_bases["contigs"] += piece.count(b"N") + piece.count(b"n")
    finally:
        for handle in handles.values():
            handle.close()
    missing = sorted(set(classes) - seen)
    if missing:
        raise ValueError(f"Sequences in the report but not in the FASTA: {missing[:5]}")
    for level in LEVELS:
        if not lengths[level]:
            raise ValueError(f"Level {level} would be empty")
        tmp[level].replace(out_dir / f"{level}.fna.gz")
    return {
        "min_gap": min_gap,
        "organelles_excluded": organelles,
        "levels": {
            level: {
                "sequences": len(lengths[level]),
                "total_length": sum(lengths[level]),
                "n50": n50(lengths[level]),
                "n_bases": n_bases[level],
            }
            for level in LEVELS
        },
    }


def make_levels(args):
    snapshot, row, run = task_context(args.work, args.task_id)
    genome = run / f"{row['accession']}.fna.gz"
    receipt = json.loads(genome.with_suffix(".gz.download.json").read_text())
    if receipt["url"] != row["fasta_url"] or receipt["sha256"] != sha256(genome):
        raise ValueError(f"Download provenance mismatch: {genome}")
    stats = build_levels(
        genome,
        read_jsonl(run / "sequence_report.jsonl"),
        run / "levels",
        snapshot["min_gap"],
    )
    stats["accession"] = row["accession"]
    stats["genome_sha256"] = receipt["sha256"]
    stats["ncbi_number_of_contigs"] = int(row["number_of_contigs"] or 0)
    save_json(run / "levels" / "levels.json", stats)
    for level, values in stats["levels"].items():
        print(f"{level}: {values['sequences']} sequences, N50 {values['n50']:,} bp")


def hobrac_command(row, fasta, output, mash_db_version, container, extra):
    for arg in extra:
        if arg.split("=", 1)[0] in FORBIDDEN_HOBRAC_ARGS:
            raise ValueError(f"HOBRAC_ARGS may not override {arg}")
    name = re.sub(r"[^A-Za-z0-9 ._-]", "_", row["organism_name"])
    command = [
        "hobrac",
        "-a",
        str(fasta),
        "-n",
        name,
        "-t",
        row["taxid"],
        "-o",
        str(output),
        "--mash-db-version",
        mash_db_version,
        "-e",
        "local",
    ]
    if container != "none":
        command.append(f"--use-{container}")
    return command + list(extra)


def host_info():
    cpu_model = ""
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                cpu_model = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    return {
        "node": platform.node(),
        "cpu_model": cpu_model,
        "cpus": len(os.sched_getaffinity(0)),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
        "slurm_mem_per_node": os.environ.get("SLURM_MEM_PER_NODE", ""),
        "slurm_partition": os.environ.get("SLURM_JOB_PARTITION", ""),
    }


def container_environment(container, work):
    env = dict(os.environ)
    if container != "none":
        key = f"{container.upper()}_BIND"
        inherited = env.get(key) or env.get(key + "PATH", "")
        env[key] = ",".join(part for part in (inherited, f"{work}:{work}") if part)
    return env


def selected_reference(output):
    path = output / "mash" / "selected_accessions.txt"
    references = [
        line.strip() for line in path.read_text().splitlines() if line.strip()
    ]
    if not references:
        raise ValueError("hobrac selected no reference")
    return references[0]


def validate_outputs(output, mash_db_version):
    version = (output / "mash" / "mash_db.msh.version").read_text().strip()
    if version != mash_db_version:
        raise ValueError(
            f"hobrac used MASH release {version}, expected {mash_db_version}"
        )
    reference = selected_reference(output)
    dotplot = output / "aln" / f"busco_{reference}" / "dotplot_busco.png"
    if not dotplot.is_file():
        raise ValueError(f"Missing BUSCO dotplot: {dotplot}")
    if not any((output / "benchmarks").glob("*.txt")):
        raise ValueError("Missing hobrac benchmarks")


def run_level(args):
    import fcntl

    snapshot, row, run = task_context(args.work, args.task_id)
    level_dir = run / args.level
    level_dir.mkdir(exist_ok=True)
    fasta = run / "levels" / f"{args.level}.fna.gz"
    with (level_dir / "task.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        state = level_dir / "status.json"
        previous = json.loads(state.read_text()) if state.exists() else {}
        if previous.get("status") == "success":
            save_json(level_dir / "completed.json", previous)
            print(f"Already complete: {row['accession']} {args.level}")
            return
        output = level_dir / "hobrac"
        if output.exists():
            # A resumed run would under-report its runtime: restart from scratch.
            shutil.rmtree(output)
        command = hobrac_command(
            row,
            fasta,
            output,
            snapshot["mash_db_version"],
            args.container,
            args.hobrac_arg,
        )
        status = {
            "status": "running",
            "attempt": previous.get("attempt", 0) + 1,
            "started_utc": now(),
            **host_info(),
        }
        save_json(state, status)
        try:
            hobrac_version = importlib.metadata.version("HoBRAC")
        except importlib.metadata.PackageNotFoundError:
            hobrac_version = "unavailable; record installation commit"
        save_json(
            level_dir / "execution.json",
            {
                "command": command,
                "command_shell": shlex.join(command),
                "hobrac_version": hobrac_version,
                "python": platform.python_version(),
                "script_sha256": sha256(__file__),
                "fasta_sha256": sha256(fasta),
            },
        )
        try:
            with (level_dir / "hobrac.log").open("a") as log:
                log.write(f"\n{now()} {shlex.join(command)}\n")
                log.flush()
                start = time.monotonic()
                subprocess.run(
                    command,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                    env=container_environment(args.container, args.work),
                )
                wall = time.monotonic() - start
            validate_outputs(output, snapshot["mash_db_version"])
            status.update(
                status="success", finished_utc=now(), hobrac_wall_seconds=round(wall, 3)
            )
            save_json(state, status)
            save_json(level_dir / "completed.json", status)
        except Exception as exc:
            status.update(status="failed", finished_utc=now(), error=str(exc))
            save_json(state, status)
            raise


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def step_group(stem):
    if stem == "busco_assembly":
        return "busco_assembly"
    if stem.startswith("busco_reference_"):
        return "busco_reference"
    if stem.startswith("aln_"):
        return "minimap2"
    if stem in ("mash", "select_references"):
        return "mash"
    if stem in ("download_db", "get_busco_dataset", "download_busco_dataset") or (
        stem.startswith("get_reference_")
    ):
        return "download"
    return "other"


def read_benchmarks(directory):
    rows = []
    for path in sorted(Path(directory).glob("*.txt")):
        records = read_tsv(path)
        if not records:
            continue
        record = records[-1]
        rows.append(
            {
                "rule": path.stem,
                "group": step_group(path.stem),
                "s": number(record.get("s")),
                "cpu_time": number(record.get("cpu_time")),
                "max_rss_mb": number(record.get("max_rss")),
                "mean_load": number(record.get("mean_load")),
            }
        )
    return rows


def busco_summary(directory):
    for path in sorted(Path(directory).rglob("short_summary*.txt")):
        text = path.read_text()
        values = {
            key: re.search(rf"\b{key}:([0-9.]+)", text)
            for key in ("C", "S", "D", "F", "M", "n")
        }
        if all(values.values()):
            return {key: float(match[1]) for key, match in values.items()}
    return {}


def mash_distance(path, reference):
    best = None
    with Path(path).open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 5 and fields[0].rpartition(":")[0] == reference:
                candidate = (float(fields[2]), float(fields[3]))
                best = candidate if best is None else min(best, candidate)
    return best


def anchor_metrics(paf, significance):
    """Share of BUSCO anchors drawn in a significant (coloured) dotplot pair."""
    counts = Counter()
    with Path(paf).open() as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) >= 6:
                counts[(fields[0], fields[5])] += 1
    significant = {
        (row["query"], row["target"])
        for row in read_tsv(significance)
        if row["significant"].strip().lower() == "true"
    }
    anchors = sum(counts.values())
    return {
        "busco_anchors": anchors,
        "anchor_sequences": len({query for query, _ in counts}),
        "significant_pairs": len(significant),
        "significant_anchor_fraction": (
            round(sum(v for k, v in counts.items() if k in significant) / anchors, 4)
            if anchors
            else ""
        ),
    }


def collect_run(work, manifest_row, level, levels_stats):
    run = work / "runs" / manifest_row["accession"]
    level_dir = run / level
    row = {
        **{
            key: manifest_row[key]
            for key in ("task_id", "accession", "organism_name", "phylum", "taxid")
        },
        "level": level,
        "status": "not_run",
    }
    if levels_stats:
        row.update(levels_stats["levels"][level])
        if level == "contigs":
            row["ncbi_number_of_contigs"] = levels_stats.get(
                "ncbi_number_of_contigs", ""
            )
    benchmarks = []
    state = level_dir / "status.json"
    if not state.exists():
        return row, benchmarks
    status = json.loads(state.read_text())
    for key in (
        "status",
        "error",
        "attempt",
        "node",
        "cpu_model",
        "cpus",
        "slurm_job_id",
        "hobrac_wall_seconds",
    ):
        row[key] = status.get(key, "")
    if status["status"] != "success":
        return row, benchmarks
    output = level_dir / "hobrac"
    try:
        reference = selected_reference(output)
        row["reference_accession"] = reference
        distance = mash_distance(output / "mash" / "mash.dist", reference)
        if distance:
            row["reference_distance"], row["reference_pvalue"] = distance
        phylum_file = output / "mash" / "mash_db.msh.phylum"
        if phylum_file.exists():
            row["hobrac_phylum"] = phylum_file.read_text().strip()
        dataset = output / "busco" / "chosen_dataset.txt"
        if dataset.exists():
            row["busco_lineage"] = dataset.read_text().split()[0]
        busco = busco_summary(output / "busco" / "busco_assembly")
        for key, field in (
            ("C", "busco_complete"),
            ("S", "busco_single"),
            ("D", "busco_duplicated"),
            ("F", "busco_fragmented"),
            ("M", "busco_missing"),
            ("n", "busco_markers"),
        ):
            row[field] = busco.get(key, "")
        aln = output / "aln" / f"busco_{reference}"
        row.update(
            anchor_metrics(aln / "aln_busco.paf", aln / "significance_results.txt")
        )
        for variant, name in (
            ("dotplot_busco", "dotplot_busco.png"),
            ("dotplot_busco_bw", "dotplot_busco_bw.png"),
        ):
            if (aln / name).is_file():
                row[variant] = str((aln / name).relative_to(work))
        benchmarks = read_benchmarks(output / "benchmarks")
        totals = defaultdict(float)
        for bench in benchmarks:
            totals[bench["group"]] += bench["s"] or 0.0
        for group in STEP_GROUPS:
            row[f"{group}_seconds"] = round(totals[group], 3)
    except (OSError, ValueError, KeyError) as exc:
        row["status"], row["error"] = "invalid_result", str(exc)
    return row, [
        {"accession": row["accession"], "level": level, **bench} for bench in benchmarks
    ]


def summarize(args):
    work = args.work
    snapshot = load_snapshot(work, verify_all=True)
    manifest = read_tsv(work / "manifest.tsv")
    report = work / "report"
    report.mkdir(exist_ok=True)
    rows, benchmarks = [], []
    for manifest_row in manifest:
        stats_path = (
            work / "runs" / manifest_row["accession"] / "levels" / "levels.json"
        )
        stats = json.loads(stats_path.read_text()) if stats_path.exists() else None
        for level in LEVELS:
            row, bench = collect_run(work, manifest_row, level, stats)
            rows.append(row)
            benchmarks.extend(bench)
    write_tsv(report / "runs.tsv", rows, RUN_FIELDS)
    write_tsv(
        report / "benchmarks.tsv",
        benchmarks,
        (
            "accession",
            "level",
            "rule",
            "group",
            "s",
            "cpu_time",
            "max_rss_mb",
            "mean_load",
        ),
    )

    selection = []
    for manifest_row in manifest:
        by_level = {
            r["level"]: r for r in rows if r["accession"] == manifest_row["accession"]
        }
        entry = {
            key: manifest_row[key] for key in ("accession", "organism_name", "phylum")
        }
        for level in LEVELS:
            entry[f"reference_{level}"] = by_level[level].get("reference_accession", "")
            entry[f"distance_{level}"] = by_level[level].get("reference_distance", "")
        references = {entry[f"reference_{level}"] for level in LEVELS}
        complete = all(by_level[level]["status"] == "success" for level in LEVELS)
        entry["same_reference"] = (
            ("yes" if len(references) == 1 else "no") if complete else "incomplete"
        )
        selection.append(entry)
    write_tsv(report / "reference_selection.tsv", selection, tuple(selection[0]))

    statuses = Counter(r["status"] for r in rows)
    save_json(
        report / "summary.json",
        {
            "created_utc": now(),
            "genomes": len(manifest),
            "runs": len(rows),
            "statuses": dict(statuses),
            "same_reference": dict(Counter(e["same_reference"] for e in selection)),
            "mash_db_version": snapshot["mash_db_version"],
            "script_sha256": sha256(__file__),
            "snapshot_sha256": sha256(work / "snapshot.json"),
        },
    )
    plot_runtime(rows, report)
    plot_signal(rows, report)
    plot_dotplot_grids(rows, work, report, args.grid_rows)
    print(json.dumps(dict(statuses), indent=2))
    print(f"Tables and plots: {report}")


LEVEL_MARKERS = {"contigs": "o", "scaffolds": "s", "chromosomes": "^"}
RUNTIME_PANELS = (
    ("End-to-end hobrac run", "hobrac_wall_seconds"),
    ("BUSCO (assembly)", "busco_assembly_seconds"),
    ("BUSCO (reference)", "busco_reference_seconds"),
    ("Minimap2 alignment", "minimap2_seconds"),
    ("MASH search", "mash_seconds"),
    ("Downloads (network)", "download_seconds"),
    ("Other steps", "other_seconds"),
)
RUNTIME_AXES = (
    ("n50", "N50 (Mb)", 1e-6, "runtime_vs_n50", "N50"),
    (
        "total_length",
        "Cumulative length (Mb)",
        1e-6,
        "runtime_vs_total_length",
        "cumulative length",
    ),
    (
        "sequences",
        "Number of sequences",
        1,
        "runtime_vs_sequences",
        "number of sequences",
    ),
)


def phylum_colors(rows):
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("tab20")
    return {p: cmap(i % 20) for i, p in enumerate(sorted({r["phylum"] for r in rows}))}


def plain_log_formatter(label_only_base=True):
    """Matplotlib's log-tick thinning, printed as 0.25 instead of 2.5e-01."""
    from matplotlib.ticker import LogFormatter

    class PlainLogFormatter(LogFormatter):
        def __call__(self, x, pos=None):
            return f"{x:.3g}" if super().__call__(x, pos) else ""

    return PlainLogFormatter(labelOnlyBase=label_only_base)


def log_axis(axis):
    axis.set_major_formatter(plain_log_formatter())
    axis.set_minor_formatter(plain_log_formatter(label_only_base=False))


def successful_genomes(rows):
    genomes = defaultdict(list)
    for row in rows:
        if row["status"] == "success":
            genomes[row["accession"]].append(row)
    return genomes


def draw_series(ax, genomes, colors, x_of, y_of):
    """One line per genome through its levels, markers by level, colours by phylum."""
    for group in genomes.values():
        points = sorted(
            (
                (x, y, r)
                for r in group
                for x, y in [(x_of(r), y_of(r))]
                if x is not None and y is not None
            ),
            key=lambda point: point[:2],
        )
        if len(points) > 1:
            ax.plot(
                [p[0] for p in points],
                [p[1] for p in points],
                color=colors[group[0]["phylum"]],
                linewidth=1,
                alpha=0.6,
            )
        for x, y, r in points:
            ax.scatter(
                x,
                y,
                marker=LEVEL_MARKERS[r["level"]],
                color=colors[r["phylum"]],
                edgecolor="black",
                linewidth=0.3,
                s=32,
                zorder=3,
            )


def draw_legends(ax, colors):
    from matplotlib.lines import Line2D

    ax.axis("off")
    phylum_legend = ax.legend(
        handles=[
            Line2D([], [], color=color, linewidth=4, label=phylum)
            for phylum, color in colors.items()
        ],
        title="Phylum (colour)",
        loc="upper left",
        ncol=2,
        fontsize=7,
        title_fontsize=8,
        labelspacing=0.3,
        handlelength=1.5,
        columnspacing=1.0,
    )
    ax.add_artist(phylum_legend)
    ax.legend(
        handles=[
            Line2D([], [], color="#444444", marker=marker, linestyle="", label=level)
            for level, marker in LEVEL_MARKERS.items()
        ],
        title="Assembly level (shape)",
        loc="lower left",
        fontsize=7,
        title_fontsize=8,
        labelspacing=0.3,
    )


def plot_runtime(rows, report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = phylum_colors(rows)
    genomes = successful_genomes(rows)
    runs = sum(len(group) for group in genomes.values())
    for key, label, scale, stem, name in RUNTIME_AXES:
        fig, axes = plt.subplots(4, 2, figsize=(8.27, 11.69))
        for ax, (title, field) in zip(axes.flat, RUNTIME_PANELS):
            draw_series(
                ax,
                genomes,
                colors,
                lambda r, key=key: positive(r, key, scale),
                lambda r, field=field: positive(r, field, 1 / 60),
            )
            ax.set_xscale("log")
            ax.set_yscale("log")
            log_axis(ax.xaxis)
            log_axis(ax.yaxis)
            ax.set_title(title, fontsize=10)
            ax.set_xlabel(label)
            ax.set_ylabel("Wall-clock time (min)")
            ax.grid(alpha=0.2, which="both")
        draw_legends(axes.flat[-1], colors)
        fig.suptitle(f"hobrac runtime per step vs {name} ({runs} runs)", fontsize=11)
        fig.tight_layout()
        for suffix in ("png", "pdf"):
            fig.savefig(report / f"{stem}.{suffix}", dpi=200)
        plt.close(fig)


def positive(row, field, scale=1.0):
    """Scaled field value, or None when missing or non-positive (log axes)."""
    value = number(row.get(field))
    return value * scale if value and value > 0 else None


def anchors_per_sequence(row):
    anchors, sequences = positive(row, "busco_anchors"), positive(row, "anchor_sequences")
    return anchors / sequences if anchors and sequences else None


def plot_signal(rows, report):
    """Share of anchors in significant pairs, by level and against contiguity."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = phylum_colors(rows)
    genomes = successful_genomes(rows)
    runs = sum(len(group) for group in genomes.values())
    fraction = lambda r: number(r.get("significant_anchor_fraction"))
    fig, axes = plt.subplots(2, 2, figsize=(8.27, 7.5))
    panels = (
        ("Assembly level", lambda r: LEVELS.index(r["level"]), None),
        ("N50 (Mb)", lambda r: positive(r, "n50", 1e-6), "log"),
        ("Anchors per anchor-bearing sequence", anchors_per_sequence, "log"),
    )
    for ax, (label, x_of, scale) in zip(axes.flat, panels):
        draw_series(ax, genomes, colors, x_of, fraction)
        if scale:
            ax.set_xscale(scale)
            log_axis(ax.xaxis)
        else:
            ax.set_xticks(range(len(LEVELS)), LEVELS)
            ax.set_xlim(-0.4, len(LEVELS) - 0.6)
        ax.set_ylim(-0.03, 1.03)
        ax.set_xlabel(label)
        ax.set_ylabel("Anchors in significant pairs (fraction)")
        ax.grid(alpha=0.2, which="both")
    draw_legends(axes.flat[-1], colors)
    fig.suptitle(f"Synteny signal vs assembly contiguity ({runs} runs)", fontsize=11)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(report / f"signal_vs_contiguity.{suffix}", dpi=200)
    plt.close(fig)


def tile_title(level, row):
    title = level
    if number(row.get("sequences")):
        title += f": {int(row['sequences']):,} seq, N50 {int(row['n50']) / 1e6:.2f} Mb"
    if row.get("reference_accession"):
        title += f"\nreference {row['reference_accession']}"
    return title


def plot_dotplot_grids(rows, work, report, rows_per_page):
    """One genome per row, levels as columns; paged PNGs plus one PDF per variant."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    by_key = {(r["accession"], r["level"]): r for r in rows}
    genomes = sorted({(r["phylum"], r["organism_name"], r["accession"]) for r in rows})
    size = rows_per_page or len(genomes)
    pages = [genomes[i : i + size] for i in range(0, len(genomes), size)]
    out_dir = report / "dotplot_grids"
    out_dir.mkdir(exist_ok=True)
    for stale in out_dir.glob("busco*_page*.png"):
        stale.unlink()
    for field, suffix in (("dotplot_busco", ""), ("dotplot_busco_bw", "_bw")):
        with PdfPages(report / f"dotplot_grid_busco{suffix}.pdf") as pdf:
            for page_number, page in enumerate(pages, 1):
                fig, axes = plt.subplots(
                    len(page),
                    len(LEVELS),
                    figsize=(4.2 * len(LEVELS) + 0.8, 4.5 * len(page)),
                    squeeze=False,
                )
                for i, (phylum, organism, accession) in enumerate(page):
                    for j, level in enumerate(LEVELS):
                        ax = axes[i][j]
                        ax.axis("off")
                        row = by_key[(accession, level)]
                        image = row.get(field)
                        if image and (work / image).is_file():
                            ax.imshow(
                                mpimg.imread(work / image), interpolation="antialiased"
                            )
                        else:
                            ax.text(
                                0.5,
                                0.5,
                                row["status"].replace("_", " "),
                                ha="center",
                                va="center",
                                fontsize=12,
                                color="#666666",
                                transform=ax.transAxes,
                            )
                        ax.set_title(tile_title(level, row), fontsize=8)
                    axes[i][0].text(
                        -0.04,
                        0.5,
                        f"{organism}\n{phylum} · {accession}",
                        transform=axes[i][0].transAxes,
                        rotation=90,
                        ha="right",
                        va="center",
                        fontsize=9,
                    )
                fig.tight_layout()
                fig.savefig(
                    out_dir / f"busco{suffix}_page{page_number:02d}.png", dpi=150
                )
                pdf.savefig(fig, dpi=150)
                plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser(
        "prepare", help="Freeze metadata and draw one genome per phylum"
    )
    prep.add_argument("--seed", type=int, required=True)
    prep.add_argument("--bioproject", default=BIOPROJECT)
    prep.add_argument("--mash-db-version", default="latest")
    prep.set_defaults(func=prepare)
    fetch = sub.add_parser(
        "download", help="Download one genome and its sequence report"
    )
    fetch.set_defaults(func=download_one)
    levels = sub.add_parser(
        "levels", help="Write the contig/scaffold/chromosome FASTAs"
    )
    levels.set_defaults(func=make_levels)
    run = sub.add_parser("run", help="Run hobrac on one genome level")
    run.add_argument("--level", choices=LEVELS, required=True)
    run.add_argument(
        "--container", choices=("none", "apptainer", "singularity"), default="none"
    )
    run.add_argument(
        "--hobrac-arg",
        action="append",
        default=[],
        help="Extra hobrac option, repeatable; use --hobrac-arg=--flag",
    )
    run.set_defaults(func=run_level)
    for command in (fetch, levels, run):
        command.add_argument("--task-id", type=int, required=True)
    check = sub.add_parser("check", help="Verify the frozen snapshot before scheduling")
    check.set_defaults(func=check_work)
    summary = sub.add_parser(
        "summarize", help="Tables, runtime plots and dotplot grids"
    )
    summary.add_argument(
        "--grid-rows", type=int, default=6, help="Genomes per grid page (0: one page)"
    )
    summary.set_defaults(func=summarize)
    for command in (prep, fetch, levels, run, check, summary):
        command.add_argument(
            "--work", type=lambda p: Path(p).expanduser().resolve(), required=True
        )
    args = parser.parse_args()
    # hobrac interpolates paths into shell commands without consistent quoting.
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
