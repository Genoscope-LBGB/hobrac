# Mash reference-availability survey

This workflow samples species across hobrac-supported eukaryotic phyla and uses Mash to find and measure their nearest eligible reference genome. The preparation step freezes both the sample manifest and reference sketches.

```text
                    download A → Mash A ─┐
Frozen manifest ─── download B → Mash B ─┼─ report
                    download C → Mash C ─┘
```

## Quick start

To run the workflow, you will need Linux, Bash, Python >=3.11, installed hobrac and host-side Mash, Snakemake 9 with its Slurm executor plugin, and matplotlib. If needed, you can install the workflow dependencies directly into your hobrac environment:

```bash
python3 -m pip install 'snakemake>=9,<10' snakemake-executor-plugin-slurm 'matplotlib>=3.6'
```

Copy **this whole directory** to your machine, making sure to include `scripts/` and `profiles/`.

### 1. Configure

```bash
cp config.example.sh config.local.sh
# Edit config.local.sh before submitting anything.
```

In `config.local.sh`, set `WORK` to an absolute shared path without spaces or shell metacharacters, choose your sample `COUNT` (which defaults to 1000), and fill in your environment and Slurm settings. `MAX_CONCURRENT` limits submitted child jobs, including pending jobs, `MASTER_TIME` sets the runtime limit for the controller job, and `TIME` limits the runtime of each Mash job.

**A different count, seed, scope, or reference snapshot requires a new `WORK` directory.**

### 2. Prepare and review

```bash
bash submit.sh prepare
```

Monitor `$WORK/logs/` and confirm that the preparation job completes successfully. Inspect `manifest.tsv`, `phylum_coverage.tsv`, and `sampling_exclusions.json` to verify the sampling results before you proceed.

### 3. Run

```bash
bash submit.sh workflow
```

By default (`CLEANUP=true`), the workflow deletes downloaded and prepared FASTA files for successful queries to save disk space. It retains all final results, mappings, logs, provenance metadata, and shared reference data. If a query fails, its FASTA files are kept so you can investigate the failure.

### 4. Inspect results

Start with these reports under `$WORK/report/`:

- `distance_by_phylum.png` / `.pdf`: Distance distributions grouped by query phylum.
- `nearest_references.tsv`: Every sampled query, including failures.
- `phylum_summary.tsv` and `summary.json`: Counts and status summaries.

## Restart or inspect a run

```bash
bash submit.sh workflow                  # Continue unfinished work
bash submit.sh summarize                 # Produce a partial report separately
bash workflow.sh config.local.sh --dry-run
```

Snakemake reuses results from completed queries automatically, even after their temporary FASTA files have been removed. If the controller times out or is killed forcefully, make sure all previous child jobs have finished or been cancelled before you restart. Never run two controllers against the same `WORK` directory. Clear a stale lock only after confirming that no controller or child job is still writing to `WORK`.

To run the workflow locally instead of through Slurm:

```bash
bash workflow.sh config.local.sh --executor local --cores 8
```

## Directory layout

| Path | Purpose |
| --- | --- |
| `config.example.sh` | Template configuration to copy and edit for your run |
| `submit.sh` | Submits preparation, the controller, or reporting jobs to Slurm |
| `workflow.sh`, `Snakefile` | Launches the Snakemake controller and defines the DAG |
| `profiles/` | Default scheduler configuration profiles |
| `scripts/` | Internal Python backend and Slurm job entry point |
| `tests/` | Offline unit tests and local workflow integration tests |

To run the test suite from this directory:

```bash
python3 -m unittest discover -s tests -t .. -v
```

These tests run entirely offline using a mock hobrac/Mash setup and a loopback FASTA server, so they do not submit real Slurm jobs or download files from NCBI. They verify step overlap, temporary file cleanup, restarts, failure handling, migration, and include a 1,000-genome DAG dry run.
