import argparse
import os
import sys


def get_args():
    parser = argparse.ArgumentParser(
        prog="HoBRAC",
        description="\n\nHomology-based reference genome acquisition and comparison",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=True,
    )

    parser.add_argument(
        "-a",
        "--assembly",
        action="store",
        dest="assembly",
        help="Path to a genome assembly fasta file",
        required=True,
        default=None,
        type=os.path.abspath,
    )
    parser.add_argument(
        "-n",
        "--name",
        action="store",
        dest="scientific_name",
        help="Scientific name of the organism given in --assembly",
        required=True,
        default=None,
    )
    parser.add_argument(
        "-t",
        "--taxid",
        action="store",
        dest="taxid",
        help="Taxid of the organism given in --assembly",
        required=True,
        default=None,
    )
    parser.add_argument(
        "--metaeuk",
        action="store_true",
        dest="metaeuk",
        help="Use metaeuk instead of miniprot",
        required=False,
        default=False,
    )

    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument(
        "-r",
        "--reference",
        action="store",
        dest="reference",
        help="Path to a fasta file to use as a reference (disables the reference searching step)",
        default=None,
        type=os.path.abspath,
    )
    optional_args.add_argument(
        "--busco-assembly",
        action="store",
        dest="busco_assembly",
        help="Path to a BUSCO result directory (or run dir or full_table.tsv) to reuse for the assembly",
        default=None,
        type=os.path.abspath,
    )
    optional_args.add_argument(
        "--busco-reference",
        action="store",
        dest="busco_reference",
        help="Path to a BUSCO result directory (or run dir or full_table.tsv) to reuse for the reference",
        default=None,
        type=os.path.abspath,
    )

    optional_args.add_argument(
        "-e",
        "--executor",
        action="store",
        dest="executor",
        help="Name of a snakemake executor plugin to execute the pipeline on a computing cluster",
        default="slurm",
    )
    optional_args.add_argument(
        "--profile",
        action="store",
        dest="profile",
        help="Path to a snakemake profile directory",
        default=None,
        type=os.path.abspath,
    )
    optional_args.add_argument(
        "--rerun-incomplete",
        action="store_true",
        dest="rerun_incomplete",
        help="Restart incomplete jobs (typically after a crash)",
    )
    optional_args.add_argument(
        "--busco-memory",
        action="store",
        dest="busco_memory",
        help="Amount of RAM in GB reserved for Busco",
        default=100,
        type=int,
    )
    optional_args.add_argument(
        "--minimap2-memory",
        action="store",
        dest="minimap2_memory",
        help="Amount of RAM in GB reserved for Minimap2",
        default=100,
        type=int,
    )
    optional_args.add_argument(
        "--busco-runtime",
        action="store",
        dest="busco_runtime",
        help="Maximum runtime in hours for Busco jobs",
        default=24,
        type=int,
    )
    optional_args.add_argument(
        "--minimap2-runtime",
        action="store",
        dest="minimap2_runtime",
        help="Maximum runtime in hours for Minimap2 jobs",
        default=12,
        type=int,
    )
    optional_args.add_argument(
        "--qos",
        action="store",
        dest="qos",
        help="SLURM QoS for all jobs",
        default=None,
    )
    optional_args.add_argument(
        "-o",
        "--output",
        action="store",
        dest="output_directory",
        help="Output directory",
        default="hobrac_analysis",
    )
    optional_args.add_argument(
        "--allow-same-taxid",
        action="store_true",
        dest="allow_same_taxid",
        help="Allows the chosen reference to be of the same taxid as the assembly",
        default=False,
        required=False,
    )
    optional_args.add_argument(
        "--allow-zero-distance",
        action="store_true",
        dest="allow_zero_distance",
        help="Allows the choice of the closest reference to pick a genome with a Mash distance of 0",
        default=False,
        required=False,
    )
    optional_args.add_argument(
        "--stop-after-mash",
        action="store_true",
        dest="stop_after_mash",
        help="Stop the pipeline after computing mash distances",
        default=False,
        required=False,
    )
    optional_args.add_argument(
        "--use-apptainer",
        action="store_true",
        dest="use_apptainer",
        help="Use apptainer to run the pipeline",
        default=False,
    )
    optional_args.add_argument(
        "--use-singularity",
        action="store_true",
        dest="use_singularity",
        help="Use singularity to run the pipeline",
        default=False,
    )
    optional_args.add_argument(
        "--use-docker",
        action="store_true",
        dest="use_docker",
        help="Use docker to run the pipeline",
        default=False,
    )

    args = parser.parse_args()

    return args
