import argparse
import os


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
        action="append",
        dest="reference",
        help=(
            "Path to a fasta file to use as a reference"
            " (disables the reference searching step)."
            " Can be specified multiple times."
        ),
        default=None,
        type=os.path.abspath,
    )
    optional_args.add_argument(
        "--busco-assembly",
        action="store",
        dest="busco_assembly",
        help=(
            "Path to a BUSCO result directory"
            " (or run dir or full_table.tsv)"
            " to reuse for the assembly"
        ),
        default=None,
        type=os.path.abspath,
    )
    optional_args.add_argument(
        "--busco-reference",
        action="store",
        dest="busco_reference",
        help=(
            "Path to a BUSCO result directory"
            " (or run dir or full_table.tsv)"
            " to reuse for the reference"
        ),
        default=None,
        type=os.path.abspath,
    )

    optional_args.add_argument(
        "-e",
        "--executor",
        action="store",
        dest="executor",
        help=(
            "Name of a snakemake executor plugin"
            " to execute the pipeline on a computing cluster"
        ),
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
        help=(
            "Allows the choice of the closest reference"
            " to pick a genome with a Mash distance of 0"
        ),
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
        "--ref-count",
        action="store",
        dest="ref_count",
        help="Number of reference genomes to select automatically",
        default=1,
        type=int,
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
        "--min-busco-genes",
        action="store",
        dest="min_busco_genes",
        help="Minimum complete BUSCO genes required per sequence for JCVI plot",
        default=30,
        type=int,
    )
    color_group = optional_args.add_mutually_exclusive_group()
    color_group.add_argument(
        "--jcvi-custom-colors",
        action="store",
        dest="jcvi_custom_colors",
        help=(
            "Path to custom color file for JCVI synteny plot"
            " (tab-separated: BUSCO_ID, COLOR, ALG_NAME;"
            " COLOR may be R,G,B, #rrggbb, or rrggbb)."
            " By default, ALG statistical testing still runs"
            " to determine significance;"
            " use --jcvi-skip-alg to disable it."
            " Genes not in the file and genes not significantly associated"
            " (unless --jcvi-skip-alg is used)"
            " will be shown in grey."
            " The ALG_NAME column is also used for"
            " rearrangement index calculation."
        ),
        default=None,
        type=os.path.abspath,
    )
    color_group.add_argument(
        "--jcvi-color-metazoan-alg",
        action="store_true",
        dest="jcvi_color_metazoan_alg",
        help=(
            "Use pre-computed 29-ALG color scheme for JCVI synteny plot."
            " The color file is automatically selected based on the BUSCO"
            " dataset. If the dataset is not part of the pre-computed set,"
            " the pipeline falls back to default coloring."
        ),
        default=False,
    )
    color_group.add_argument(
        "--jcvi-color-bilaterian-alg",
        action="store_true",
        dest="jcvi_color_bilaterian_alg",
        help=(
            "Use pre-computed 24-bilaterian-ALG color scheme for JCVI"
            " synteny plot. The color file is automatically selected"
            " based on the BUSCO dataset. If the dataset is not part"
            " of the pre-computed set, the pipeline falls back to"
            " default coloring."
        ),
        default=False,
    )
    optional_args.add_argument(
        "--jcvi-names",
        action="store",
        dest="jcvi_names",
        help="Comma-separated custom names for JCVI tracks. If one name is provided, "
        "it applies only to the assembly. If multiple names are provided, the count "
        "must equal 1 (assembly) + number of references, in order.",
        default="",
    )
    optional_args.add_argument(
        "--jcvi-hide-non-significant",
        action="store_true",
        dest="hide_non_significant",
        help="Hide links between chromosome pairs without significant associations "
        "in the JCVI karyotype plot. This produces a cleaner plot showing only "
        "ALG-related synteny.",
        default=False,
    )
    optional_args.add_argument(
        "--jcvi-skip-alg",
        action="store_true",
        dest="skip_alg",
        help="Skip ALG statistical testing when using custom colors. All genes in the "
        "color file get their custom color; unlisted genes are grey. Only effective "
        "with --jcvi-custom-colors.",
        default=False,
    )
    optional_args.add_argument(
        "--jcvi-pvalue",
        action="store",
        dest="jcvi_pvalue",
        help="Base significance threshold (alpha) for ALG detection via Fisher's "
        "exact test. Bonferroni correction is applied on top of this value.",
        default=0.01,
        type=float,
    )
    optional_args.add_argument(
        "--jcvi-min-chain-genes",
        action="store",
        dest="jcvi_min_chain_genes",
        help=(
            "Minimum BUSCO genes a chromosome chain must be supported by to"
            " appear in the final output. Chains with fewer genes are dropped"
            " before sub-chain pruning, so sub-chains hidden inside a dropped"
            " long chain can re-emerge."
        ),
        default=5,
        type=int,
    )

    args = parser.parse_args()

    return args
