import argparse
import os
import sys


def get_args():
    parser = argparse.ArgumentParser(
        prog="HoBRAC",
        description="\n\nHomology-based reference genome acquisition and comparison",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
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

    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument(
        "-e",
        "--executor",
        action="store",
        dest="executor",
        help="Name of a snakemake executor plugin to execute the pipeline on a computing cluster",
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
        help="Allow the chosen reference to be of the same taxid as the assembly",
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
    optional_args.add_argument(
        "--container-version",
        action="store_true",
        dest="container_version",
        help="Container version to use",
        default="0.1",
    )

    args = parser.parse_args()

    return args
