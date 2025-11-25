#!/usr/bin/env python3
import glob
import os
import shutil
import subprocess
import sys
from hobrac.command_line import get_args


thisdir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
snakefile_path = os.path.join(thisdir, "workflow", "Snakefile")


def check_dependencies(require_busco: bool = True, require_reference_search: bool = True):
    deps = ["taxonkit", "mash"]
    if require_reference_search:
        deps.extend(["find_reference_genomes", "datasets"])
    if require_busco:
        deps.append("busco")

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


def _busco_dir_contains_full_table(dir_path: str) -> bool:
    pattern = os.path.join(dir_path, "run*", "full_table.tsv")
    return len(glob.glob(pattern)) > 0


def normalize_busco_dir(path: str) -> str:
    """Normalize a user-provided BUSCO path to the directory that contains run*/full_table.tsv.
    Accepts:
      - a full_table.tsv path
      - a run_* directory path
      - a directory that contains run*/full_table.tsv
    Returns an absolute path to the directory that contains run*/full_table.tsv, or exits on error.
    """
    p = os.path.abspath(path)
    if os.path.isfile(p) and os.path.basename(p) == "full_table.tsv":
        p = os.path.dirname(p)
    # If the path itself is a run_* directory, move up one level
    if os.path.isdir(p) and os.path.basename(p).startswith("run_"):
        parent = os.path.dirname(p)
        if _busco_dir_contains_full_table(parent):
            return parent
    # Otherwise, verify the directory contains run*/full_table.tsv
    if os.path.isdir(p) and _busco_dir_contains_full_table(p):
        return p

    print(f"Provided BUSCO path does not contain run*/full_table.tsv: {path}", file=sys.stderr)
    sys.exit(1)


def link_busco_dir(src_dir: str, dest_dir: str):
    """Create a symlink from dest_dir to src_dir. Fails if a non-symlink exists at dest_dir."""
    create_dir(os.path.dirname(dest_dir))
    if os.path.lexists(dest_dir):
        if os.path.islink(dest_dir):
            os.unlink(dest_dir)
        else:
            print(
                f"Destination {dest_dir} already exists and is not a symlink. Remove it or choose a new output directory.",
                file=sys.stderr,
            )
            sys.exit(1)
    os.symlink(src_dir, dest_dir, target_is_directory=True)


def get_base_snakemake_args(args) -> str:
    cmd = "snakemake --latency-wait 70 --jobs 100 -p --skip-script-cleanup "

    # Fallback to greedy scheduler if CBC solver is not available to avoid PuLP errors
    if shutil.which("cbc") is None:
        cmd += "--scheduler greedy "

    if args.rerun_incomplete:
        cmd += "--rerun-incomplete "

    if args.executor:
        cmd += f"--executor {args.executor} --cores 4000 "

    if args.executor == "slurm":
        cmd += "--slurm-keep-successful-logs "

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
    cmd += f"stop_after_mash={args.stop_after_mash} "

    if args.metaeuk:
        cmd += "busco_method=metaeuk "
    else:
        cmd += "busco_method=miniprot "

    cmd += f"minimap2_memory={args.minimap2_memory * 1000} "
    cmd += f"busco_memory={args.busco_memory * 1000} "

    if getattr(args, "busco_assembly_override_path", None):
        cmd += f"busco_assembly_override='{args.busco_assembly_override_path}' "
    if getattr(args, "busco_reference_override_path", None):
        cmd += f"busco_reference_override='{args.busco_reference_override_path}' "

    return cmd


def main():
    args = get_args()

    create_dir(args.output_directory)
    os.chdir(args.output_directory)

    # Handle precomputed BUSCO results: normalize and symlink
    if getattr(args, "busco_assembly", None):
        norm = normalize_busco_dir(args.busco_assembly)
        link_busco_dir(norm, os.path.join("busco", "busco_assembly"))
        args.busco_assembly_override_path = os.path.abspath(os.path.join("busco", "busco_assembly"))
    else:
        args.busco_assembly_override_path = None

    if getattr(args, "busco_reference", None):
        norm = normalize_busco_dir(args.busco_reference)
        link_busco_dir(norm, os.path.join("busco", "busco_reference"))
        args.busco_reference_override_path = os.path.abspath(os.path.join("busco", "busco_reference"))
    else:
        args.busco_reference_override_path = None

    # Dependencies: require busco only if at least one side still needs to run
    require_busco = not (args.busco_assembly_override_path and args.busco_reference_override_path)
    check_dependencies(
        require_busco=require_busco,
        require_reference_search=not args.stop_after_mash,
    )

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
