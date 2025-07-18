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
    url: str | None = None


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="precompute_mash_refseq",
        description="\n\nAutomatically download reference genomes and compute their MASH sketches",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
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

    get_eukaryote_list(args.output_dir)
    extract_chromosome_complete(args.output_dir)
    download_taxdump(args.output_dir)
    extract_phylum(args.output_dir)
    phylums = collect_phylums(args.output_dir)

    download_dir = os.path.join(args.output_dir, "downloads")
    if not os.path.exists(download_dir):
        os.mkdir(download_dir)
        
    phylums = get_ncbi_genome_ftp_url(phylums, download_dir)    
    download_genomes(phylums, download_dir)
    decompress(download_dir)

    mash_dir = os.path.join(args.output_dir, "mash")
    mash_dir_tmp = os.path.join(mash_dir, "tmp")
    if not os.path.exists(mash_dir_tmp):
        os.makedirs(mash_dir_tmp, exist_ok=True)
    run_mash(download_dir, mash_dir_tmp)
    paste_mash(mash_dir_tmp, mash_dir)


def get_eukaryote_list(output_dir):
    print("Getting genome list...", flush=True, file=sys.stderr)
    
    ftp_url = "ftp.ncbi.nih.gov"
    ftp = ftplib.FTP(ftp_url)
    ftp.login()

    url = "genomes/GENOME_REPORTS/eukaryotes.txt"
    print(f"Dowloading ftp://{ftp_url}/{url}", file=sys.stderr)

    with open(f"{output_dir}/eukaryotes.txt", "wb") as f:
        ftp.retrbinary(f"RETR {url}", f.write)
    ftp.quit()


def extract_chromosome_complete(output_dir):
    print("Extracting complete genomes...", flush=True, file=sys.stderr)
    
    with open(f"{output_dir}/eukaryotes.txt") as inf, open(f"{output_dir}/eukaryotes_chr_complete.txt", "w") as out:
        for line in inf:
            if "Chromosome" in line or "Complete Genome" in line:
                line = line.split("\t")
                print(f"{line[1]}\t{line[8]}", file=out)
                
                
def download_taxdump(output_dir):
    print("Downloading taxdump...", flush=True, file=sys.stderr)
    
    taxdump = os.path.join(output_dir, "taxdump")
    os.makedirs(taxdump, exist_ok=True)
    
    taxdump_gz = os.path.join(taxdump, "taxdump.tar.gz")
    os.system(f"wget -c https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O {taxdump_gz}")
    os.system(f"cd {taxdump} && tar -xzf {taxdump_gz} && rm {taxdump_gz}")
    

def extract_phylum(output_dir):
    print("Extracting phylums...", flush=True, file=sys.stderr)
    
    genome_list = os.path.join(output_dir, "eukaryotes_chr_complete.txt")
    taxdump = os.path.join(output_dir, "taxdump")
    accession_list = os.path.join(output_dir, "final_list.txt")
    
    os.environ["TAXONKIT_DB"] = taxdump
    os.system(f"cat {genome_list} | taxonkit reformat -I 1 --format '{{p}}' -r 'no_returned_phylum' " \
        f"| grep -v 'no_returned_phylum' > {accession_list}")
    
    shutil.rmtree(taxdump)
    
    
def collect_phylums(output_dir) -> defaultdict[str, List[Genome]]:
    accession_list = os.path.join(output_dir, "final_list.txt")
    
    phylums: defaultdict[str, List[Genome]] = defaultdict(list)
    with open(accession_list) as inf:
        for line in inf:
            line = line.rstrip("\n").split("\t")
            accession: str = line[1]
            phylum: str = line[2].replace(" ", "_")
            phylums[phylum].append(Genome(accession, None))
            
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


def download_genomes(phylums: defaultdict[str, List[Genome]], output_dir):
    already_downloaded_path = os.path.join(output_dir, "already_downloaded.txt")

    for phylum in phylums:
        phylums_download_dir = os.path.join(output_dir, phylum)
        os.makedirs(phylums_download_dir, exist_ok=True)
        
        for genome in phylums[phylum]:
            if genome.url is None:
                continue
            
            accession = genome.accession
            url = genome.url
            compressed_name = os.path.join(phylums_download_dir, f"{accession}.fna.gz")

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


def decompress(output_dir):
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


def run_mash(input_dir, output_dir):
    pattern = r"(?P<phylum>[^/]+)/(?P<accession>GC[AF]_\d{9}\.\d+)"
    
    for f in glob.glob(f"{input_dir}/*/*.fna"):
        match = re.search(pattern, f)
        if match is None:
            print(f"Could not find a GC[AF] match for {f}", flush=True, file=sys.stderr)
            continue
        
        phylum = match.group("phylum")
        identifier = match.group("accession")
        output_path = os.path.join(output_dir, phylum, os.path.basename(f))
        os.makedirs(output_path, exist_ok=True)

        subprocess.run(["mash", "sketch", "-s", "10000", "-I", identifier, "-o", output_path, f], check=True)
        os.remove(f)


def paste_mash(input_dir, output_dir):
    for f in glob.glob(f"{output_dir}/*.fofn"):
        os.remove(f)
    
    pattern = r"(?P<phylum>[^/]+)/(?P<accession>GC[AF]_\d{9}\.\d+)"
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

    shutil.rmtree(input_dir)
