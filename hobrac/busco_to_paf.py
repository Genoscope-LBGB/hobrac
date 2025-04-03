#!/usr/bin/env python3
import argparse
import os
import sys
import uuid
from Bio import SeqIO

def read_busco_tsv(file_path):
    busco_data = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if parts[1] == "Complete":
                busco_id = parts[0]
                busco_data[busco_id] = {
                    'chr': parts[2],
                    'start': int(parts[3]),
                    'end': int(parts[4])
                }
    return busco_data

def calculate_fasta_lengths(file_path):
    lengths = {}
    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, "fasta"):
            lengths[record.id] = len(record.seq)
    return lengths

def write_idx_file(file_path, keyword, lengths):
    with open(file_path, 'w') as file:
        file.write(f"{keyword}\n")
        for seq_id, length in lengths.items():
            file.write(f"{seq_id}\t{length}\n")

def generate_paf(busco1_data, busco2_data, len_query, len_target, output_file):
    with open(output_file, 'w') as outfile:
        for busco_id, data in busco1_data.items():
            if busco_id in busco2_data:
                target = busco2_data[busco_id]['chr']
                outfile.write(f"{data['chr']}\t{len_query[data['chr']]}\t{data['start']}\t{data['end']}\t+\t{target}\t{len_target[target]}\t{busco2_data[busco_id]['start']}\t{busco2_data[busco_id]['end']}\t1000\t1000\t60\n")

def run(busco_query, busco_ref, query_fasta, ref_fasta, output_dir=None):
    # Generate a unique output directory name if not provided
    if output_dir is None:
        output_dir = f"alignment_busco_{uuid.uuid4()}"

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    busco_query_data = read_busco_tsv(busco_query)
    busco_ref_data = read_busco_tsv(busco_ref)
    len_query = calculate_fasta_lengths(query_fasta)
    len_target = calculate_fasta_lengths(ref_fasta)

    # Write PAF file
    paf_file = os.path.join(output_dir, "aln_busco.paf")
    generate_paf(busco_query_data, busco_ref_data, len_query, len_target, paf_file)

    # Write query_assembly.idx file
    query_idx_file = os.path.join(output_dir, "query_assembly.idx")
    write_idx_file(query_idx_file, "Assembly", len_query)

    # Write target_reference.idx file
    target_idx_file = os.path.join(output_dir, "target_reference.idx")
    write_idx_file(target_idx_file, "Reference", len_target)


def main():
    parser = argparse.ArgumentParser(description="Generate PAF and index files from BUSCO and FASTA files.")
    parser.add_argument('--busco_query', required=True, help="Path to the query BUSCO TSV file.")
    parser.add_argument('--busco_ref', required=True, help="Path to the reference BUSCO TSV file.")
    parser.add_argument('--query', required=True, help="Path to the query FASTA file.")
    parser.add_argument('--ref', required=True, help="Path to the reference FASTA file.")
    parser.add_argument('--out', required=False, help="Output directory. If not specified, a unique directory name will be generated.")

    args = parser.parse_args()

    run(args.busco_query, args.busco_ref, args.query, args.ref, args.out)
