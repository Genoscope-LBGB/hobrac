#!/usr/bin/env bash
# Internal Slurm entry point. Submit through submit.sh, not directly.
set -euo pipefail
# Slurm spools this script, so submit.sh must pass the real experiment root.
STAGE=${1:?stage required}
CONFIG=${2:?absolute config path required}
ROOT=${3:?absolute experiment root required}
# The controller launcher loads its own configuration and environment once.
if [[ $STAGE == workflow ]]; then
    exec bash "$ROOT/workflow.sh" "$CONFIG"
fi
source "$CONFIG"
if [[ -n ${ENV_SETUP:-} ]]; then source "$ENV_SETUP"; fi
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1

case "$STAGE" in
    prepare)
        python3 "$ROOT/scripts/fragmentation.py" prepare --work "$WORK" \
            --seed "$SEED" --bioproject "$BIOPROJECT" --mash-db-version "$MASH_DB_VERSION"
        ;;
    summarize)
        python3 "$ROOT/scripts/fragmentation.py" summarize --work "$WORK" \
            --grid-rows "${GRID_ROWS:-6}"
        ;;
    *) echo "Unknown stage: $STAGE" >&2; exit 2 ;;
esac
