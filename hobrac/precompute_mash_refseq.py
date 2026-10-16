import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass

import requests

NCBI_URL = "https://ftp.ncbi.nlm.nih.gov"
# GenBank groups holding eukaryotes; their assembly_summary.txt files give the
# FTP directory of every current assembly in a few seconds of download.
GENBANK_GROUPS = (
    "fungi",
    "invertebrate",
    "plant",
    "protozoa",
    "vertebrate_mammalian",
    "vertebrate_other",
)
TIMEOUT = 60
DOWNLOAD_ATTEMPTS = 3

_local = threading.local()


@dataclass
class Genome:
    accession: str
    taxid: str
    phylum: str
    url: str | None = None


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="precompute_mash_refseq",
        description=(
            "\n\nAutomatically download reference genomes"
            " and compute their MASH sketches. Rerunning with the same output"
            " directory only downloads genomes that are not sketched yet."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-e",
        "--eukaryote-list",
        action="store",
        dest="euk_list_file",
        help=(
            "File containing the list of"
            " chromosome/complete eukaryote genomes"
            " at NCBI. TSV format:"
            " GCA_accession\\tAssembly_name\\tTaxid"
        ),
        required=True,
        type=os.path.abspath,
    )
    parser.add_argument(
        "-o",
        action="store",
        dest="output_dir",
        help=(
            "Output directory. Per-genome sketches are kept in sketches/ and"
            " the per-phylum databases are written to mash/"
        ),
        required=True,
        type=os.path.abspath,
    )
    parser.add_argument(
        "-t",
        "--threads",
        action="store",
        dest="threads",
        help="Number of genomes downloaded and sketched concurrently",
        default=4,
        type=int,
    )

    return parser.parse_args()


def main():
    args = get_args()
    os.makedirs(args.output_dir, exist_ok=True)

    download_taxdump(args.output_dir)
    extract_phylum(args.output_dir, args.euk_list_file)
    genomes = read_genomes(args.output_dir)

    sketch_dir = os.path.join(args.output_dir, "sketches")
    download_dir = os.path.join(args.output_dir, "downloads")
    os.makedirs(sketch_dir, exist_ok=True)
    os.makedirs(download_dir, exist_ok=True)

    pending = [
        genome
        for genome in genomes
        if not os.path.exists(sketch_path(sketch_dir, genome))
    ]
    log(f"{len(genomes) - len(pending)} genomes already sketched, {len(pending)} to do")

    failures = resolve_urls(pending)
    failures += sketch_genomes(
        [genome for genome in pending if genome.url is not None],
        sketch_dir,
        download_dir,
        args.threads,
    )
    write_failures(args.output_dir, failures)

    paste_sketches(genomes, sketch_dir, os.path.join(args.output_dir, "mash"))


def log(message: str):
    print(message, flush=True, file=sys.stderr)


def session() -> requests.Session:
    """One HTTP session per thread, reusing connections to NCBI."""
    if not hasattr(_local, "session"):
        _local.session = requests.Session()
    return _local.session


def download_taxdump(output_dir):
    log("Downloading taxdump...")

    taxdump = os.path.join(output_dir, "taxdump")
    os.makedirs(taxdump, exist_ok=True)

    taxdump_gz = os.path.join(taxdump, "taxdump.tar.gz")
    os.system(
        f"wget -c https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O {taxdump_gz}"
    )
    os.system(f"cd {taxdump} && tar -xzf {taxdump_gz} && rm {taxdump_gz}")


def extract_phylum(output_dir, euk_list_file):
    log("Extracting phylums...")

    genome_list = euk_list_file
    taxdump = os.path.join(output_dir, "taxdump")
    accession_list = os.path.join(output_dir, "final_list.txt")

    os.environ["TAXONKIT_DB"] = taxdump
    os.system(
        f"cat {genome_list}"
        f" | taxonkit reformat -I 3"
        f" --format '{{p}}' -r 'no_returned_phylum'"
        f" > {accession_list}"
    )

    shutil.rmtree(taxdump)


def read_genomes(output_dir) -> list[Genome]:
    accession_list = os.path.join(output_dir, "final_list.txt")

    genomes = []
    with open(accession_list) as inf:
        for line in inf:
            fields = line.rstrip("\n").split("\t")
            genomes.append(Genome(fields[0], fields[2], fields[3].replace(" ", "_")))

    return genomes


def sketch_path(sketch_dir: str, genome: Genome) -> str:
    return os.path.join(sketch_dir, f"{genome.accession}_{genome.taxid}.msh")


def resolve_urls(genomes: list[Genome]) -> list[tuple[str, str]]:
    """Set each genome's FASTA URL; return (accession, reason) for failures."""
    ftp_paths = fetch_ftp_paths()

    failures = []
    for genome in genomes:
        directory = ftp_paths.get(genome.accession)
        if directory is None:
            # Not in the current GenBank eukaryote summaries, e.g. a GCF or a
            # replaced version: list its directory instead.
            try:
                directory = find_assembly_directory(genome.accession)
            except (requests.RequestException, ValueError) as e:
                failures.append((genome.accession, f"URL lookup failed: {e}"))
                continue
            if directory is None:
                failures.append((genome.accession, "assembly directory not found"))
                continue
        genome.url = f"{directory}/{directory.rsplit('/', 1)[1]}_genomic.fna.gz"

    log(f"Resolved {len(genomes) - len(failures)}/{len(genomes)} genome URLs")
    return failures


def fetch_ftp_paths() -> dict[str, str]:
    """Map accession to HTTPS assembly directory from the GenBank summaries."""
    paths = {}
    for group in GENBANK_GROUPS:
        url = f"{NCBI_URL}/genomes/genbank/{group}/assembly_summary.txt"
        log(f"Downloading {url}")
        response = session().get(url, timeout=TIMEOUT)
        response.raise_for_status()

        columns = None
        for line in response.text.splitlines():
            if line.startswith("#"):
                header = line.lstrip("#").strip().split("\t")
                if "ftp_path" in header:
                    columns = header
                continue
            if columns is None:
                raise ValueError(f"No header with an ftp_path column in {url}")
            row = dict(zip(columns, line.split("\t")))
            if row.get("ftp_path", "na") != "na":
                ftp_path = row["ftp_path"].replace("ftp://", "https://")
                paths[row["assembly_accession"]] = ftp_path.rstrip("/")
    return paths


def find_assembly_directory(accession: str) -> str | None:
    match = re.fullmatch(r"(GC[AF])_(\d{3})(\d{3})(\d{3})\.\d+", accession)
    if not match:
        raise ValueError(f"invalid accession format: {accession}")

    parent = f"{NCBI_URL}/genomes/all/{'/'.join(match.groups())}"
    response = session().get(f"{parent}/", timeout=TIMEOUT)
    if response.status_code == 404:
        return None
    response.raise_for_status()

    folder = re.search(rf'href="({re.escape(accession)}_[^"/]+)/?"', response.text)
    return f"{parent}/{folder.group(1)}" if folder else None


def sketch_genomes(
    genomes: list[Genome], sketch_dir: str, download_dir: str, threads: int
) -> list[tuple[str, str]]:
    """Download and sketch genomes concurrently; return failures."""
    failures = []
    pool = ThreadPoolExecutor(max_workers=threads)
    try:
        futures = {
            pool.submit(sketch_genome, genome, sketch_dir, download_dir): genome
            for genome in genomes
        }
        for done, future in enumerate(as_completed(futures), start=1):
            genome = futures[future]
            try:
                future.result()
                status = "sketched"
            except subprocess.CalledProcessError as e:
                failures.append((genome.accession, f"mash failed: {e.stderr.strip()}"))
                status = "mash failed"
            except (requests.RequestException, OSError) as e:
                failures.append((genome.accession, f"download failed: {e}"))
                status = "download failed"
            log(f"[{done}/{len(genomes)}] {genome.accession}: {status}")
    finally:
        # On Ctrl-C, do not keep working through the whole queue.
        pool.shutdown(cancel_futures=True)
    return failures


def sketch_genome(genome: Genome, sketch_dir: str, download_dir: str):
    """Download one compressed genome, sketch it and keep only the sketch."""
    assert genome.url is not None
    compressed = os.path.join(download_dir, f"{genome.accession}.fna.gz")
    prefix = os.path.join(download_dir, genome.accession)
    try:
        download(genome.url, compressed)
        subprocess.run(
            [
                "mash",
                "sketch",
                "-s",
                "10000",
                "-I",
                f"{genome.accession}:{genome.taxid}",
                "-o",
                prefix,
                compressed,
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        # The sketch only appears once complete, so its presence marks the
        # genome as done for later runs.
        os.replace(f"{prefix}.msh", sketch_path(sketch_dir, genome))
    finally:
        for leftover in (compressed, f"{prefix}.msh"):
            if os.path.exists(leftover):
                os.remove(leftover)


def download(url: str, path: str):
    for attempt in range(1, DOWNLOAD_ATTEMPTS + 1):
        try:
            with session().get(url, stream=True, timeout=TIMEOUT) as response:
                response.raise_for_status()
                expected = response.headers.get("Content-Length")
                with open(path, "wb") as out:
                    for chunk in response.iter_content(chunk_size=1 << 20):
                        out.write(chunk)
                    written = out.tell()
            # mash would silently sketch a truncated gzip.
            if expected is not None and written != int(expected):
                raise OSError(f"truncated download: {written}/{expected} bytes")
            return
        except (requests.RequestException, OSError) as e:
            client_error = (
                isinstance(e, requests.HTTPError)
                and e.response is not None
                and e.response.status_code < 500
            )
            if client_error or attempt == DOWNLOAD_ATTEMPTS:
                raise
            log(f"Retrying {url} after error: {e}")
            time.sleep(10 * attempt)


def write_failures(output_dir: str, failures: list[tuple[str, str]]):
    """Failures are retried on the next run; this file only reports them."""
    failed_path = os.path.join(output_dir, "failed.tsv")
    with open(failed_path, "w") as out:
        for accession, reason in failures:
            reason = " ".join(reason.split())
            print(f"{accession}\t{reason}", file=out)
    if failures:
        log(f"{len(failures)} genomes failed, see {failed_path}")


def paste_sketches(genomes: list[Genome], sketch_dir: str, mash_dir: str):
    """Rebuild one database per phylum from the sketches of the listed genomes.

    Sketches of genomes absent from the current list (e.g. replaced assembly
    versions) stay in sketch_dir but are left out of the databases.
    """
    sketches_by_phylum = defaultdict(list)
    for genome in genomes:
        sketch = sketch_path(sketch_dir, genome)
        if os.path.exists(sketch):
            sketches_by_phylum[genome.phylum].append(sketch)

    os.makedirs(mash_dir, exist_ok=True)
    for stale in os.listdir(mash_dir):
        if stale.endswith(".msh"):
            os.remove(os.path.join(mash_dir, stale))

    for phylum, sketches in sorted(sketches_by_phylum.items()):
        log(f"Pasting {len(sketches)} sketches into {phylum}.msh")
        with tempfile.NamedTemporaryFile(
            "w", dir=mash_dir, suffix=".fofn", delete=False
        ) as fofn:
            fofn.write("\n".join(sketches) + "\n")
        try:
            subprocess.run(
                [
                    "mash",
                    "paste",
                    "-l",
                    os.path.join(mash_dir, phylum),
                    fofn.name,
                ],
                check=True,
            )
        finally:
            os.remove(fofn.name)
