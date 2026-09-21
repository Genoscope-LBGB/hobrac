#!/usr/bin/env bash
# Usage: bash submit.sh prepare|array|retry|summarize [config.local.sh]
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
STAGE=${1:?Usage: bash submit.sh prepare|array|retry|summarize [config.local.sh]}
CONFIG=$(realpath "${2:-$SCRIPT_DIR/config.local.sh}")
source "$CONFIG"
if [[ -n ${ENV_SETUP:-} ]]; then source "$ENV_SETUP"; fi
if [[ ! $WORK =~ ^/[A-Za-z0-9_./-]+$ || $WORK == /path/to/* ]]; then
    echo "Set WORK to a real absolute shared path without spaces or shell metacharacters." >&2
    exit 2
fi
mkdir -p "$WORK/logs"
COMMON=(--parsable --chdir="$SCRIPT_DIR" --cpus-per-task="$CPUS"
        --output="$WORK/logs/%x-%A_%a.out")
[[ -z ${PARTITION:-} ]] || COMMON+=(--partition="$PARTITION")
[[ -z ${ACCOUNT:-} ]] || COMMON+=(--account="$ACCOUNT")
[[ -z ${QOS:-} ]] || COMMON+=(--qos="$QOS")
COMMON+=("${EXTRA_SBATCH[@]}")

submit_summary() {
    sbatch "${COMMON[@]}" --job-name=hobrac-survey-summary \
        --mem="$SUMMARY_MEMORY" --time="$SUMMARY_TIME" "$@" \
        "$SCRIPT_DIR/job.sh" summarize "$CONFIG" "$SCRIPT_DIR"
}

case "$STAGE" in
    prepare)
        sbatch "${COMMON[@]}" --job-name=hobrac-survey-prepare \
            --mem="$PREP_MEMORY" --time="$PREP_TIME" \
            "$SCRIPT_DIR/job.sh" prepare "$CONFIG" "$SCRIPT_DIR"
        echo "After success, inspect manifest.tsv and phylum_coverage.tsv, then submit the array."
        ;;
    array|retry)
        N=$(python3 - "$SCRIPT_DIR" "$WORK" <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from survey import load_snapshot, read_tsv
work = Path(sys.argv[2])
load_snapshot(work)
print(len(read_tsv(work / 'manifest.tsv')))
PY
)
        [[ $N -gt 0 && $MAX_CONCURRENT =~ ^[1-9][0-9]*$ ]]
        if [[ $STAGE == array ]]; then
            TASKS="0-$((N - 1))"
        else
            TASKS=$(<"$WORK/report/retry_task_ids.txt")
            if [[ -z $TASKS ]]; then echo "No retry tasks in the latest report."; exit 0; fi
            [[ $TASKS =~ ^[0-9]+(,[0-9]+)*$ ]] || { echo "Invalid retry task list" >&2; exit 2; }
        fi
        JOB=$(sbatch "${COMMON[@]}" --job-name=hobrac-survey-mash \
            --mem="$MEMORY" --time="$TIME" --array="${TASKS}%${MAX_CONCURRENT}" \
            "$SCRIPT_DIR/job.sh" run "$CONFIG" "$SCRIPT_DIR")
        # --parsable can return jobid;cluster on federated installations.
        JOB=${JOB%%;*}
        echo "Mash array: $JOB"
        echo "Summary job (runs even if some queries fail):"
        submit_summary --dependency="afterany:$JOB"
        ;;
    summarize) submit_summary ;;
    *) echo "Unknown stage: $STAGE" >&2; exit 2 ;;
esac
