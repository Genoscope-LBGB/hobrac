#!/usr/bin/env bash
# Run the Snakemake controller in an existing allocation (or locally).
# Usage: bash workflow.sh [config.local.sh] [extra Snakemake arguments ...]
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
CONFIG=$(realpath "${1:-$SCRIPT_DIR/config.local.sh}")
if [[ $# -gt 0 ]]; then shift; fi
source "$CONFIG"
if [[ -n ${ENV_SETUP:-} ]]; then source "$ENV_SETUP"; fi
if [[ ! $WORK =~ ^/[A-Za-z0-9_./-]+$ || $WORK == /path/to/* ]]; then
    echo "Set WORK to a real absolute shared path without spaces or shell metacharacters." >&2
    exit 2
fi
[[ $MAX_CONCURRENT =~ ^[1-9][0-9]*$ ]] || { echo "Invalid MAX_CONCURRENT" >&2; exit 2; }
export TAXONKIT_DB="$WORK/metadata/taxonomy"
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1
PROFILE=${PROFILE:-$SCRIPT_DIR/profiles/slurm}
# Fail closed before scheduling if a snapshot or successful result was modified.
python3 "$SCRIPT_DIR/scripts/survey.py" check --work "$WORK"
# A per-launch config file preserves required defaults when the user supplies
# --config overrides. Child jobs can read it from the shared WORK directory.
mkdir -p "$WORK/.snakemake"
WORKFLOW_CONFIG="$WORK/.snakemake/launch-config-$$.json"
python3 - "$WORKFLOW_CONFIG" \
    work="$WORK" container="$CONTAINER" cleanup="${CLEANUP:-true}" \
    cpus="$CPUS" memory="$MEMORY" time="$TIME" \
    download_memory="${DOWNLOAD_MEMORY:-2G}" download_time="${DOWNLOAD_TIME:-01:00:00}" \
    account="$ACCOUNT" partition="$PARTITION" qos="$QOS" \
    -- "${EXTRA_SBATCH[@]}" <<'PY'
import json
import shlex
import sys
from pathlib import Path
separator = sys.argv.index("--", 2)
settings = dict(arg.split("=", 1) for arg in sys.argv[2:separator])
# The Slurm plugin appends this resource to its shell submission command.
settings["slurm_extra"] = shlex.join(sys.argv[separator + 1:])
Path(sys.argv[1]).write_text(json.dumps(settings))
PY
# Native Snakemake CLI options supplied last can override profile settings.
exec snakemake --snakefile "$SCRIPT_DIR/Snakefile" --directory "$WORK" \
    --profile "$PROFILE" --jobs "$MAX_CONCURRENT" \
    --configfile "$WORKFLOW_CONFIG" "$@"
