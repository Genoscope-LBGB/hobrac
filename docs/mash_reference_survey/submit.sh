#!/usr/bin/env bash
# Usage: bash submit.sh prepare|workflow|retry|summarize [config.local.sh]
# "array" remains an alias for "workflow"; no large Slurm array is submitted.
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
STAGE=${1:?Usage: bash submit.sh prepare|workflow|retry|summarize [config.local.sh]}
CONFIG=$(realpath "${2:-$SCRIPT_DIR/config.local.sh}")
source "$CONFIG"
if [[ -n ${ENV_SETUP:-} ]]; then source "$ENV_SETUP"; fi
if [[ ! $WORK =~ ^/[A-Za-z0-9_./-]+$ || $WORK == /path/to/* ]]; then
    echo "Set WORK to a real absolute shared path without spaces or shell metacharacters." >&2
    exit 2
fi
mkdir -p "$WORK/logs"
COMMON=(--parsable --chdir="$SCRIPT_DIR" --cpus-per-task=1
        --output="$WORK/logs/%x-%A_%a.out")
[[ -z ${PARTITION:-} ]] || COMMON+=(--partition="$PARTITION")
[[ -z ${ACCOUNT:-} ]] || COMMON+=(--account="$ACCOUNT")
[[ -z ${QOS:-} ]] || COMMON+=(--qos="$QOS")
COMMON+=("${EXTRA_SBATCH[@]}")

submit_summary() {
    sbatch "${COMMON[@]}" --job-name=hobrac-survey-summary \
        --mem="$SUMMARY_MEMORY" --time="$SUMMARY_TIME" "$@" \
        "$SCRIPT_DIR/scripts/job.sh" summarize "$CONFIG" "$SCRIPT_DIR"
}

case "$STAGE" in
    prepare)
        sbatch "${COMMON[@]}" --job-name=hobrac-survey-prepare \
            --mem="$PREP_MEMORY" --time="$PREP_TIME" \
            "$SCRIPT_DIR/scripts/job.sh" prepare "$CONFIG" "$SCRIPT_DIR"
        echo "After success, inspect manifest.tsv and phylum_coverage.tsv, then submit the workflow."
        ;;
    workflow|array|retry)
        [[ -s $WORK/snapshot.json && -s $WORK/manifest.tsv ]] || {
            echo "Run preparation successfully before submitting the workflow." >&2; exit 2;
        }
        [[ $MAX_CONCURRENT =~ ^[1-9][0-9]*$ ]] || { echo "Invalid MAX_CONCURRENT" >&2; exit 2; }
        # Only the controller is submitted here. Snakemake fills a bounded job pool.
        sbatch "${COMMON[@]}" --job-name=hobrac-survey-workflow \
            --mem="$SUMMARY_MEMORY" --time="${MASTER_TIME:-24:00:00}" \
            "$SCRIPT_DIR/scripts/job.sh" workflow "$CONFIG" "$SCRIPT_DIR"
        ;;
    summarize) submit_summary ;;
    *) echo "Unknown stage: $STAGE" >&2; exit 2 ;;
esac
