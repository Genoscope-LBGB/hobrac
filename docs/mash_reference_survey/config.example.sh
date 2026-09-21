# Copy to config.local.sh and edit
WORK=/path/to/shared/scratch/hobrac_mash_survey

# Optional script that activates the EXISTING hobrac installation and tools.
# Example contents: source /path/to/miniforge/etc/profile.d/conda.sh; conda activate hobrac
ENV_SETUP=""

COUNT=100
SEED=20260717
# All supported eukaryotic phyla by default
ANCESTOR_TAXID=2759
# none, apptainer, or singularity
# hobrac and host-side mash must be installed
CONTAINER=none

# Optional Slurm options (leave empty to use your cluster's defaults).
PARTITION=""
ACCOUNT=""
QOS=""
# Additional sbatch options, e.g. EXTRA_SBATCH=(--constraint=...)
EXTRA_SBATCH=()
MAX_CONCURRENT=4
CPUS=2
MEMORY=24G
TIME=04:00:00
PREP_MEMORY=16G
PREP_TIME=04:00:00
SUMMARY_MEMORY=16G
SUMMARY_TIME=01:00:00

# For container runs set TAXONKIT_DB to an existing database so the launcher
# accepts it. job.sh overrides it with the survey's frozen taxonomy for run jobs.
# Native taxonkit uses the same frozen database through TAXONKIT_DB.
