#!/usr/bin/env bash
# Usage: bash submit.sh prepare|workflow|summarize [config.local.sh]
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
STAGE=${1:?Usage: bash submit.sh prepare|workflow|summarize [config.local.sh]}
CONFIG=$(realpath "${2:-$SCRIPT_DIR/config.local.sh}")
source "$CONFIG"
if [[ -n ${ENV_SETUP:-} ]]; then source "$ENV_SETUP"; fi
if [[ ! $WORK =~ ^/[A-Za-z0-9_./-]+$ || $WORK == /path/to/* ]]; then
    echo "Set WORK to a real absolute shared path without spaces or shell metacharacters." >&2
    exit 2
fi
mkdir -p "$WORK/logs"
COMMON=(--parsable --chdir="$SCRIPT_DIR" --cpus-per-task=1
        --output="$WORK/logs/%x-%j.out")
[[ -z ${PARTITION:-} ]] || COMMON+=(--partition="$PARTITION")
[[ -z ${ACCOUNT:-} ]] || COMMON+=(--account="$ACCOUNT")
[[ -z ${QOS:-} ]] || COMMON+=(--qos="$QOS")
COMMON+=("${EXTRA_SBATCH[@]}")

case "$STAGE" in
    prepare)
        sbatch "${COMMON[@]}" --job-name=hobrac-fragmentation-prepare \
            --mem="$PREP_MEMORY" --time="$PREP_TIME" \
            "$SCRIPT_DIR/scripts/job.sh" prepare "$CONFIG" "$SCRIPT_DIR"
        echo "After success, inspect manifest.tsv and phylum_coverage.tsv, then submit the workflow."
        ;;
    workflow)
        [[ -s $WORK/snapshot.json && -s $WORK/manifest.tsv ]] || {
            echo "Run preparation successfully before submitting the workflow." >&2; exit 2;
        }
        [[ $MAX_CONCURRENT =~ ^[1-9][0-9]*$ ]] || { echo "Invalid MAX_CONCURRENT" >&2; exit 2; }
        # Only the controller is submitted here. Snakemake fills a bounded job pool.
        sbatch "${COMMON[@]}" --job-name=hobrac-fragmentation-workflow \
            --mem="$SUMMARY_MEMORY" --time="$MASTER_TIME" \
            "$SCRIPT_DIR/scripts/job.sh" workflow "$CONFIG" "$SCRIPT_DIR"
        ;;
    summarize)
        sbatch "${COMMON[@]}" --job-name=hobrac-fragmentation-summary \
            --mem="$SUMMARY_MEMORY" --time="$SUMMARY_TIME" \
            "$SCRIPT_DIR/scripts/job.sh" summarize "$CONFIG" "$SCRIPT_DIR"
        ;;
    *) echo "Unknown stage: $STAGE" >&2; exit 2 ;;
esac
