#!/usr/bin/env bash
# Internal Slurm entry point. Submit through submit.sh, not directly.
set -euo pipefail
# Slurm spools this script, so submit.sh must pass the real survey root.
STAGE=${1:?stage required}
CONFIG=${2:?absolute config path required}
SURVEY_ROOT=${3:?absolute survey root required}
# The controller launcher loads its own configuration and environment once.
if [[ $STAGE == workflow ]]; then
    exec bash "$SURVEY_ROOT/workflow.sh" "$CONFIG"
fi
source "$CONFIG"
if [[ -n ${ENV_SETUP:-} ]]; then source "$ENV_SETUP"; fi
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1

case "$STAGE" in
    prepare)
        python3 "$SURVEY_ROOT/scripts/survey.py" prepare --work "$WORK" \
            --count "$COUNT" --seed "$SEED" --ancestor-taxid "$ANCESTOR_TAXID"
        ;;
    summarize)
        python3 "$SURVEY_ROOT/scripts/survey.py" summarize --work "$WORK"
        ;;
    *) echo "Unknown stage: $STAGE" >&2; exit 2 ;;
esac
