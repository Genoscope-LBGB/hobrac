# Copy to config.local.sh and edit
WORK=/path/to/shared/scratch/hobrac_fragmentation

# Optional script that activates the EXISTING hobrac installation and tools.
# Example contents: source /path/to/miniforge/etc/profile.d/conda.sh; conda activate hobrac
ENV_SETUP=""

SEED=20261002
# ATLASea umbrella BioProject
BIOPROJECT=PRJEB64126
# "latest" is resolved once during preparation and frozen in snapshot.json;
# every hobrac run then uses that dated release. YYYYMMDD pins it explicitly.
MASH_DB_VERSION=latest
# none, apptainer, or singularity
CONTAINER=none

# Optional Slurm options (leave empty to use your cluster's defaults).
PARTITION=""
ACCOUNT=""
QOS=""
# Additional sbatch options for every job. Runtimes are only comparable on one
# node type, e.g. EXTRA_SBATCH=(--constraint=...)
EXTRA_SBATCH=()
# Maximum submitted child jobs (downloads, level FASTAs and hobrac runs share it).
MAX_CONCURRENT=6
# The controller must outlive the slowest hobrac run plus queueing.
MASTER_TIME=7-00:00:00
# Optional custom Snakemake profile; empty uses profiles/slurm beside the scripts.
PROFILE=""
# Delete downloaded and level FASTAs once their hobrac runs succeed.
CLEANUP=true
DOWNLOAD_MEMORY=2G
DOWNLOAD_TIME=02:00:00
LEVELS_MEMORY=8G
LEVELS_TIME=02:00:00
# One Slurm allocation per hobrac run; hobrac executes all of its steps inside
# it (-e local). Both BUSCO jobs request 12 threads each.
CPUS=24
MEMORY=160G
TIME=3-00:00:00
# Extra hobrac options for every run, e.g. HOBRAC_ARGS=(--metaeuk)
HOBRAC_ARGS=()
# Genomes per dotplot-grid page; 0 puts every genome on one page.
GRID_ROWS=6
PREP_MEMORY=8G
PREP_TIME=02:00:00
SUMMARY_MEMORY=16G
SUMMARY_TIME=02:00:00

# For container runs, workflow.sh points TAXONKIT_DB at the frozen taxonomy.
# Optional: export NCBI_API_KEY to raise the NCBI Datasets rate limit.
