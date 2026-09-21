#!/usr/bin/env bash
# Internal Slurm entry point. Submit through submit.sh, not directly.
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# Slurm copies this script to its spool directory, so submit.sh passes the real path.
STAGE=${1:?stage required}
CONFIG=${2:?absolute config path required}
SCRIPT_DIR=${3:-$SCRIPT_DIR}
source "$CONFIG"
if [[ -n ${ENV_SETUP:-} ]]; then source "$ENV_SETUP"; fi
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1

case "$STAGE" in
    prepare)
        python3 "$SCRIPT_DIR/survey.py" prepare --work "$WORK" \
            --count "$COUNT" --seed "$SEED" --ancestor-taxid "$ANCESTOR_TAXID"
        ;;
    run)
        export TAXONKIT_DB="$WORK/metadata/taxonomy"
        python3 "$SCRIPT_DIR/survey.py" run --work "$WORK" \
            --task-id "${SLURM_ARRAY_TASK_ID:?array task required}" --container "$CONTAINER"
        ;;
    summarize)
        python3 "$SCRIPT_DIR/survey.py" summarize --work "$WORK"
        ;;
    *) echo "Unknown stage: $STAGE" >&2; exit 2 ;;
esac
