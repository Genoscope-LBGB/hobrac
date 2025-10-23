import argparse
import ftplib
import glob
import gzip
import os
import re
import requests
import shutil
import subprocess
import sys
import time

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List


@dataclass
class Genome:
    accession: str
    taxid: str
    url: str | None = None


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="precompute_mash_refseq",
        description="\n\nAutomatically download reference genomes and compute their MASH sketches",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
    )

    parser.add_argument(
        "-e",
        "--eukaryote-list",
        action="store",
        dest="euk_list_file",
        help="File containing the list of chromosome/complete eukaryote genomes at NCBI. TSV format: GCA_accession\\tAssembly_name\\tTaxid",
        required=True,
        default=None,
        type=os.path.abspath,
    )
    parser.add_argument(
        "-o",
        action="store",
        dest="output_dir",
        help="Output directory to store temporary genomes and mash sketches",
        required=True,
        default=None,
        type=os.path.abspath,
    )

    args = parser.parse_args()

    return args


def main():
    args = get_args()

    download_taxdump(args.output_dir)
    extract_phylum(args.output_dir, args.euk_list_file)
    phylums = collect_phylums(args.output_dir)

    download_dir = os.path.join(args.output_dir, "downloads")
    if not os.path.exists(download_dir):
        os.mkdir(download_dir)

    mash_dir = os.path.join(args.output_dir, "mash")
    mash_dir_tmp = os.path.join(mash_dir, "tmp")
    if not os.path.exists(mash_dir_tmp):
        os.makedirs(mash_dir_tmp, exist_ok=True)

    phylums = get_ncbi_genome_ftp_url(phylums, download_dir)

    # Download, decompress, and process each genome sequentially to minimize disk usage
    download_and_process_genomes(phylums, download_dir, mash_dir_tmp)

    # Paste all mash sketches together
    paste_mash(mash_dir_tmp, mash_dir)


def download_taxdump(output_dir):
    print("Downloading taxdump...", flush=True, file=sys.stderr)

    taxdump = os.path.join(output_dir, "taxdump")
    os.makedirs(taxdump, exist_ok=True)

    taxdump_gz = os.path.join(taxdump, "taxdump.tar.gz")
    os.system(f"wget -c https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O {taxdump_gz}")
    os.system(f"cd {taxdump} && tar -xzf {taxdump_gz} && rm {taxdump_gz}")


def extract_phylum(output_dir, euk_list_file):
    print("Extracting phylums...", flush=True, file=sys.stderr)

    genome_list = euk_list_file
    taxdump = os.path.join(output_dir, "taxdump")
    accession_list = os.path.join(output_dir, "final_list.txt")

    os.environ["TAXONKIT_DB"] = taxdump
    os.system(f"cat {genome_list} | taxonkit reformat -I 3 --format '{{p}}' -r 'no_returned_phylum' " f" > {accession_list}")

    shutil.rmtree(taxdump)


def collect_phylums(output_dir) -> defaultdict[str, List[Genome]]:
    accession_list = os.path.join(output_dir, "final_list.txt")

    phylums: defaultdict[str, List[Genome]] = defaultdict(list)
    with open(accession_list) as inf:
        for line in inf:
            line = line.rstrip("\n").split("\t")
            accession: str = line[0]
            taxid: str = line[2]
            phylum: str = line[3].replace(" ", "_")
            phylums[phylum].append(Genome(accession, taxid, None))

    return phylums


def get_ncbi_genome_ftp_url(phylums: defaultdict[str, List[Genome]], output_dir):
    ftp: ftplib.FTP = ftplib.FTP()

    try:
        base_https_url = "https://ftp.ncbi.nih.gov"

        ftp = ftplib.FTP("ftp.ncbi.nih.gov", timeout=30)
        ftp.login()

        already_downloaded_path = os.path.join(output_dir, "already_downloaded.txt")

        downloaded_list = set()
        if os.path.exists(already_downloaded_path):
            with open(already_downloaded_path) as inf:
                for line in inf:
                    line = line.rstrip("\n")
                    downloaded_list.add(line)

        with open(already_downloaded_path, "a") as current_list:
            for phylum in phylums:
                for genome in phylums[phylum]:
                    gca_accession = genome.accession

                    if gca_accession in downloaded_list:
                        continue

                    match = re.match(r"GC[AF]_(\d{3})(\d{3})(\d{3})", gca_accession)
                    if not match:
                        print(f"Invalid GCA accession format: {gca_accession}")
                        print(gca_accession, file=current_list)
                        continue

                    partial_path = f"/genomes/all/GCA/{match.group(1)}/{match.group(2)}/{match.group(3)}"
                    try:
                        ftp.cwd(partial_path)
                    except:
                        print(gca_accession, file=current_list)
                        continue

                    folders = []
                    ftp.retrlines("LIST", lambda line: folders.append(line))
                    full_folder_name = None
                    for folder_line in folders:
                        parts = folder_line.split()
                        if parts and parts[-1].startswith(gca_accession):
                            full_folder_name = parts[-1]
                            break

                    if not full_folder_name:
                        print(f"Could not find folder for {gca_accession}, url: {partial_path}")
                        print(gca_accession, file=current_list)
                        continue

                    ftp.cwd(full_folder_name)
                    file_name = f"{full_folder_name}_genomic.fna.gz"
                    genome.url = f"{base_https_url}{partial_path}/{full_folder_name}/{file_name}"

                    time.sleep(0.5)

    except Exception as e:
        print(e)
        try:
            ftp.quit()
        except:
            pass
        return get_ncbi_genome_ftp_url(phylums, output_dir)

    return phylums


def download_and_process_genomes(phylums: defaultdict[str, List[Genome]], download_dir: str, mash_output_dir: str):
    """Download genomes and immediately process them: decompress → mash sketch → delete."""
    already_downloaded_path = os.path.join(download_dir, "already_downloaded.txt")

    for phylum in phylums:
        phylums_download_dir = os.path.join(download_dir, phylum)
        os.makedirs(phylums_download_dir, exist_ok=True)

        for genome in phylums[phylum]:
            if genome.url is None:
                continue

            accession = genome.accession
            taxid = genome.taxid
            url = genome.url
            compressed_name = os.path.join(phylums_download_dir, f"{accession}_{taxid}.fna.gz")

            # Download the compressed genome
            try:
                print(f"Downloading {accession} from {url}", flush=True, file=sys.stderr)
                response = requests.get(url, stream=True)
                response.raise_for_status()
            except Exception as e:
                print(f"Error downloading {accession}: {e}", flush=True, file=sys.stderr)
                continue

            with open(compressed_name, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)

            # Mark as downloaded
            with open(already_downloaded_path, "a") as current_list:
                print(accession, file=current_list)

            # Immediately decompress the genome
            try:
                decompressed_path = decompress_single_file(compressed_name)
            except Exception as e:
                print(f"Error decompressing {accession}: {e}, skipping mash processing", flush=True, file=sys.stderr)
                continue

            # Immediately run mash sketch on the decompressed genome
            try:
                run_mash_single_file(decompressed_path, phylum, accession, taxid, mash_output_dir)
            except Exception as e:
                print(f"Error running mash on {accession}: {e}", flush=True, file=sys.stderr)
                # Clean up decompressed file if mash fails
                if os.path.exists(decompressed_path):
                    os.remove(decompressed_path)
                continue


def download_genomes(phylums: defaultdict[str, List[Genome]], output_dir):
    """Download genomes only (legacy function, not used in optimized workflow)."""
    already_downloaded_path = os.path.join(output_dir, "already_downloaded.txt")

    for phylum in phylums:
        phylums_download_dir = os.path.join(output_dir, phylum)
        os.makedirs(phylums_download_dir, exist_ok=True)

        for genome in phylums[phylum]:
            if genome.url is None:
                continue

            accession = genome.accession
            taxid = genome.taxid
            url = genome.url
            compressed_name = os.path.join(phylums_download_dir, f"{accession}_{taxid}.fna.gz")

            try:
                response = requests.get(url, stream=True)
                response.raise_for_status()
            except Exception as e:
                print(e)
                continue

            with open(compressed_name, "wb") as f, open(already_downloaded_path, "a") as current_list:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)

                print(accession, file=current_list)


def decompress_single_file(compressed_path: str) -> str:
    """Decompress a single .gz file and return the path to the decompressed file."""
    print(f"Decompressing {compressed_path}", flush=True, file=sys.stderr)

    decompressed_path = compressed_path.replace(".gz", "")
    chunk_size = 10 * 1024 * 1024  # 10 MB

    try:
        with gzip.open(compressed_path, "rb") as f_in:
            with open(decompressed_path, "wb") as f_out:
                while True:
                    chunk = f_in.read(chunk_size)
                    if not chunk:
                        break
                    f_out.write(chunk)
    except Exception as e:
        print(f"Error decompressing {compressed_path}: {e}", flush=True, file=sys.stderr)
        raise

    # Remove compressed file after successful decompression
    os.remove(compressed_path)
    return decompressed_path


def decompress(output_dir):
    """Decompress all .gz files in the output directory (legacy function, not used in optimized workflow)."""
    for f in glob.glob(f"{output_dir}/*/*.gz"):
        print(f"Decompressing {f}")

        chunk_size = 10 * 1024 * 1024  # 10 MB
        try:
            with gzip.open(f, "rb") as f_in:
                with open(f.replace(".gz", ""), "wb") as f_out:
                    while True:
                        chunk = f_in.read(chunk_size)
                        if not chunk:
                            break
                        f_out.write(chunk)

        except Exception as e:
            print(e)
            pass

        os.remove(f)


def run_mash_single_file(fna_path: str, phylum: str, accession: str, taxid: str, output_dir: str) -> str:
    """Run mash sketch on a single decompressed genome file and delete it afterward."""
    identifier = f"{accession}:{taxid}"
    output_path = os.path.join(output_dir, phylum, os.path.basename(fna_path))
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    print(f"Running mash sketch on {fna_path}", flush=True, file=sys.stderr)
    subprocess.run(["mash", "sketch", "-s", "10000", "-I", identifier, "-o", output_path, fna_path], check=True)

    # Remove decompressed file after mash processing
    os.remove(fna_path)
    print(f"Deleted {fna_path} after mash processing", flush=True, file=sys.stderr)

    return f"{output_path}.msh"


def run_mash(input_dir, output_dir):
    """Run mash sketch on all .fna files (legacy function, not used in optimized workflow)."""
    pattern = r".*/(?P<phylum>[^/]+)/(?P<accession>GC[AF]_\d{9}\.\d+)_(?P<taxid>\d+)\.fna$"

    for f in glob.glob(f"{input_dir}/*/*.fna"):
        match = re.search(pattern, f)
        if match is None:
            print(f"Could not find a GC[AF] match for {f}", flush=True, file=sys.stderr)
            continue

        accession = match.group("accession")
        phylum = match.group("phylum")
        taxid = match.group("taxid")
        identifier = f"{accession}:{taxid}"
        output_path = os.path.join(output_dir, phylum, os.path.basename(f))
        os.makedirs(output_path, exist_ok=True)

        subprocess.run(["mash", "sketch", "-s", "10000", "-I", identifier, "-o", output_path, f], check=True)
        os.remove(f)


def paste_mash(input_dir, output_dir):
    for f in glob.glob(f"{output_dir}/*.fofn"):
        os.remove(f)

    pattern = r".*/(?P<phylum>[^/]+)/(?P<accession>GC[AF]_\d{9}\.\d+)_(?P<taxid>\d+)\.fna.msh$"
    phylums = set()
    for f in glob.glob(f"{input_dir}/*/*.msh"):
        match = re.search(pattern, f)
        if match is None:
            print(f"Could not find a GC[AF] match for {f}", flush=True, file=sys.stderr)
            continue
        phylum = match.group("phylum")
        phylums.add(phylum)

        mash_fofn = os.path.join(output_dir, f"{phylum}.fofn")
        with open(mash_fofn, "a") as out:
            print(f, file=out)

    for phylum in phylums:
        mash_fofn = os.path.join(output_dir, f"{phylum}.fofn")
        subprocess.run(["mash", "paste", f"{input_dir}/{phylum}.msh", "-l", mash_fofn], check=True)
        os.remove(mash_fofn)

        global_sketch = os.path.join(output_dir, f"{phylum}.msh")
        if os.path.exists(global_sketch):
            new_global_sketch = os.path.join(output_dir, "new_references")
            subprocess.run(["mash", "paste", new_global_sketch, global_sketch, f"{input_dir}/{phylum}.msh"], check=True)
            os.rename(f"{new_global_sketch}.msh", global_sketch)
        else:
            os.rename(f"{input_dir}/{phylum}.msh", global_sketch)

    # shutil.rmtree(input_dir)
