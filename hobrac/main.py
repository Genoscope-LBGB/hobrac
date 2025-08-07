#!/usr/bin/env python3
import glob
import os
import shutil
import subprocess
import sys
from hobrac.command_line import get_args


thisdir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
snakefile_path = os.path.join(thisdir, "workflow", "Snakefile")


def check_dependencies():
    deps = ["find_reference_genomes", "datasets", "taxonkit", "busco", "mash"]

    for dep in deps:
        if not shutil.which(dep):
            print(f"{dep} not found, exiting.", file=sys.stderr)
            exit(1)


def create_dir(path: str):
    try:
        os.mkdir(path)
    except FileExistsError:
        pass
    except FileNotFoundError:
        print(f"Path to {path} does not exist.")
        exit(-1)
    except PermissionError:
        print(f"Unsufficient permissions to write output directory {path}")
        exit(-1)


def get_base_snakemake_args(args) -> str:
    cmd = "snakemake --latency-wait 30 --jobs 100 -p "
    
    if args.rerun_incomplete:
        cmd += "--rerun-incomplete "

    if args.executor:
        cmd += f"--executor {args.executor} --cores 4000 "

    if args.use_apptainer:
        cmd += "--use-apptainer "
    elif args.use_singularity:
        cmd += "--use-singularity "
    elif args.use_docker:
        cmd += "--use_docker "

    return cmd


def generate_snakemake_command(args) -> str:
    cmd = get_base_snakemake_args(args)

    cmd += f"--snakefile {snakefile_path} "

    cmd += "--config "
    cmd += f"assembly={args.assembly} "
    cmd += f"scientific_name='{args.scientific_name}' "
    cmd += f"taxid={args.taxid} "
    cmd += f"allow_same_taxid={args.allow_same_taxid} "
    cmd += f"allow_zero_distance={args.allow_zero_distance} "
    
    if args.miniprot:
        cmd += "busco_method=miniprot "
    else:
        cmd += "busco_method=metaeuk "
    
    cmd += f"minimap2_memory={args.minimap2_memory} "
    cmd += f"busco_memory={args.busco_memory} "

    cmd += f"container_version='{args.container_version}' "

    return cmd


def main():
    check_dependencies()

    args = get_args()

    create_dir(args.output_directory)
    os.chdir(args.output_directory)
    
    if args.reference:
        skip_reference_search(args.reference)

    cmd = generate_snakemake_command(args)
    print(f"\n{cmd}\n")

    process = subprocess.Popen(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
    process.wait()

    sys.exit(process.returncode)


def skip_reference_search(reference: str):
    create_dir("mash")
    with open("mash/closest_reference.txt", "w") as out:
        print(os.path.basename(reference), file=out, end="")
        
    create_dir("reference")
    with open("reference/reference.txt", "w") as out:
        print(reference, file=out, end="")