# Mash reference-availability survey

This workflow samples species across the eukaryotic phyla that hobrac supports and uses Mash to find and measure their nearest eligible reference genome. The preparation step freezes both the sample manifest and the reference sketches.

```text
                    download A → Mash A ─┐
Frozen manifest ─── download B → Mash B ─┼─ report
                    download C → Mash C ─┘
```

## Results and interpretation

The archived run evaluated **1,000 distinct species from 41 eukaryotic phyla** using Mash 2.3 and inputs frozen on 22 September 2026. Queries were selected from 9,942 eligible species with chromosome-level or complete assemblies, with one assembly per species. Sampling was balanced across phyla where possible: all eligible species were included in 30 phyla, while 61–62 species were sampled from each of 11 larger pools.

![Nearest eligible reference Mash distances by query phylum, colored by taxonomic relationship and shaped by Mash signal](report/outputs/distance_by_phylum.png)

[Download the PDF](report/outputs/distance_by_phylum.pdf).

Each point represents one query's nearest eligible reference from another species. A lower Mash distance indicates greater sequence similarity. Color shows the taxonomic relationship, and shape distinguishes P ≤ 0.01, P > 0.01, and no shared hashes. Phylum labels show queries with an eligible reference / sampled. The legend percentages use the 994 queries with an eligible reference.

### Reference availability and taxonomic proximity

All 1,000 queries completed, and there were no failed or incomplete jobs. **994 had an eligible reference** after same-taxid references, zero-distance matches, the query assembly, and other same-species references were excluded. References with unresolved species assignments were also excluded.

| Relationship to the selected reference | Genomes | Percentage of the 994 with a reference |
| --- | ---: | ---: |
| Same genus | 503 | 50.6% |
| Same family, different genus | 212 | 21.3% |
| More distant than family | 279 | 28.1% |

The remaining six queries were from **Acanthocephala, Chytridiomycota, Discosea, Endomyxa, Glomeromycota, and Phoronida**. Each of these phyla had only one eligible query species, and the same-taxid or zero-distance rules excluded all of its reference entries. The `0/1` labels therefore indicate no eligible alternative for the sampled species, not necessarily a database with only one assembly: Endomyxa had two reference entries and Glomeromycota had six. Under this survey's database and rules, these queries could not be evaluated against a different-species reference. This is not a general limitation of hobrac when another suitable reference can be supplied.

Among the 994 selected query–reference pairs:

- **990** had shared hashes with **P ≤ 0.01**.
- **3** had shared hashes with **P > 0.01**.
- **1** had **no shared hashes**. Its distance of 1 is a no-signal placeholder, not a measured evolutionary distance.

All three sampled Fornicata had weak or undetectable signal: *Giardia duodenalis* and *Giardia muris* each shared only 1 of 10,000 sketch hashes with their selected reference, while *Spironucleus salmonicida* shared none. The other weak match was *Phymatolithon calcareum* (Rhodophyta), with 2 of 10,000 shared hashes.

A small Mash P-value indicates that the sketch overlap is unlikely under Mash's random-sequence model. It does **not** establish that the reference is close enough for reliable whole-genome alignment or structural comparison, and these are pairwise Mash P-values, not a test of differences between phyla.

### Do larger reference databases give more same-genus matches?

The following examples compare reference database size with the selected reference's relationship and distance. Reference entries are the number of comparison rows in the phylum database **before exclusions**, not the number of unique species or the size of the eligible query pool. Same-genus percentages use the queries with an eligible reference in that phylum as the denominator.

| Phylum | Reference entries | Queries with a reference | Same-genus matches | Median Mash distance |
| --- | ---: | ---: | ---: | ---: |
| Streptophyta | 4,362 | 61 | 36 (59.0%) | 0.067 |
| Chordata | 4,221 | 61 | 29 (47.5%) | 0.033 |
| Arthropoda | 3,894 | 62 | 35 (56.5%) | 0.110 |
| Ascomycota | 1,768 | 61 | 38 (62.3%) | 0.174 |
| Apicomplexa | 150 | 45 | 40 (88.9%) | 0.129 |
| Euglenozoa | 85 | 44 | 34 (77.3%) | 0.027 |
| Bacillariophyta | 19 | 25 | 8 (32.0%) | 0.291 |
| Rhodophyta | 17 | 19 | 4 (21.1%) | 0.296 |

Some sparsely represented groups have relatively distant matches, but **database size alone does not consistently predict the same-genus match rate**. Euglenozoa, for example, has fewer reference entries than Chordata but a higher same-genus rate. Taxonomic composition matters: 25 of the 44 Euglenozoa queries are *Leishmania* and another six are *Trypanosoma*. Sampling one assembly per species avoids repeated species, but it does not balance genera or families.

Euglenozoa has the lowest median distance (0.027), but the distribution is broad: its middle 50% spans **0.005–0.201**, compared with **0.018–0.078** in Chordata. Rhodophyta has a much higher median (0.296), although that median is not the highest across all phyla. Both differences between groups and variation within groups matter when you assess reference availability. Results based on only a few species should not be generalized to an entire phylum.

### Conclusion

Finding a same-genus reference for **50.6%** of queries with an eligible reference is encouraging evidence of taxonomic proximity. Neither this statistic nor low Mash distances, however, demonstrate that hobrac selects the **best** reference for downstream analysis.

Here, the selection criterion is explicitly **minimum Mash distance**, a proxy for sequence similarity. The genome that is most similar in sequence is not necessarily the most similar in chromosome structure or gene order. Establishing which reference is best for a particular application would require comparing candidate references using downstream alignment, synteny, or structural-conservation results. This survey stops after Mash and does not perform that validation.

Reference availability depends on the genomes present in the database. Future sequencing and database updates may improve coverage, but this single snapshot cannot measure or guarantee that improvement. The appropriate conclusion is therefore **broad availability of detectable sequence matches, with uneven reference proximity**, rather than demonstrated optimal reference selection.

### Archived inputs and outputs

The files supporting these conclusions are organized as follows:

```text
report/
├── inputs/
│   ├── eligible_species.tsv
│   ├── manifest.tsv
│   ├── phylum_coverage.tsv
│   ├── sampling_exclusions.json
│   └── snapshot.json
└── outputs/
    ├── nearest_references.tsv
    ├── phylum_summary.tsv
    ├── summary.json
    ├── distance_by_phylum.png
    ├── distance_by_phylum.pdf
    ├── retry_task_ids.txt
    └── reset_required_task_ids.txt
```

See the [sample manifest](report/inputs/manifest.tsv), [sampling coverage](report/inputs/phylum_coverage.tsv), and [per-query results](report/outputs/nearest_references.tsv) for the underlying data. The input files and original provenance records are unchanged. The figures were regenerated from the per-query results with the revised colors, shapes, and legend counts. The script hashes in the provenance records identify the original run, not the later plotting changes.

This is a **report archive**, not a complete runnable work directory. It does not include the frozen reference sketches, taxonomy files, raw Mash comparisons, or execution logs. Paths inside `snapshot.json` refer to the original work-directory layout. The input/output split here does not change the workflow layout described below: preparation still writes inputs under `$WORK/`, and reporting writes outputs under `$WORK/report/`.

# Running the workflow

## Quick start

To run the workflow, you need Linux, Bash, Python >=3.11, installed hobrac and host-side Mash, Snakemake 9 with its Slurm executor plugin, and matplotlib. If needed, you can install the workflow dependencies directly into your hobrac environment:

```bash
python3 -m pip install 'snakemake>=9,<10' snakemake-executor-plugin-slurm 'matplotlib>=3.6'
```

Copy **this whole directory** to your machine, and make sure the copy includes `scripts/` and `profiles/`.

### 1. Configure

```bash
cp config.example.sh config.local.sh
# Edit config.local.sh before submitting anything.
```

In `config.local.sh`, set `WORK` to an absolute shared path without spaces or shell metacharacters, choose your sample `COUNT` (which defaults to 1000), and fill in your environment and Slurm settings. `MAX_CONCURRENT` limits submitted child jobs, including pending jobs. `MASTER_TIME` sets the runtime limit for the controller job, and `TIME` limits the runtime of each Mash job.

**A different count, seed, scope, or reference snapshot requires a new `WORK` directory.**

### 2. Prepare and review

```bash
bash submit.sh prepare
```

Monitor `$WORK/logs/` and confirm that the preparation job completes successfully. Before you proceed, inspect `manifest.tsv`, `phylum_coverage.tsv`, and `sampling_exclusions.json` and verify the sampling results.

### 3. Run

```bash
bash submit.sh workflow
```

By default (`CLEANUP=true`), the workflow deletes downloaded and prepared FASTA files for successful queries so that it can save disk space. It retains all final results, mappings, logs, provenance metadata, and shared reference data. If a query fails, its FASTA files are kept so that you can investigate the failure.

### 4. Inspect results

Start with these reports under `$WORK/report/`:

- `distance_by_phylum.png` / `.pdf`: Distance distributions grouped by query phylum.
- `nearest_references.tsv`: Every sampled query, including failures.
- `phylum_summary.tsv` and `summary.json`: Counts and status summaries.

In the distance plot, point color shows the nearest reference's taxonomic relationship: **blue** for the same genus (or closer), **orange** for the same family but a different genus, and **reddish purple** for a relationship more distant than family. The color legend includes genome counts and percentages for each category, calculated over all queries with an eligible reference, including weak matches and queries with no shared hashes. Its title states this denominator, and failed, incomplete, and no-eligible-reference queries are excluded. Point shape shows the Mash signal: **circles** for shared hashes with P ≤ 0.01, **triangles** for shared hashes with P > 0.01, and **crosses** for no shared hashes. Crosses retain the relationship color, and their distance of 1 is a placeholder for no detected overlap, not a measured evolutionary distance. Boxplots include only queries with shared hashes. These display choices do not change reference selection or the reported statistics.

## Restart or inspect a run

```bash
bash submit.sh workflow                  # Continue unfinished work
bash submit.sh summarize                 # Produce a partial report separately
bash workflow.sh config.local.sh --dry-run
```

Snakemake reuses results from completed queries automatically, even after their temporary FASTA files have been removed. If the controller times out or is killed forcefully, make sure all previous child jobs have finished or been cancelled before you restart. Never run two controllers against the same `WORK` directory. Clear a stale lock only after you confirm that no controller or child job is still writing to `WORK`.

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
| `report/inputs/` | Archived sampling inputs and snapshot provenance |
| `report/outputs/` | Archived result tables, figures, and status summaries |

To run the test suite from this directory:

```bash
python3 -m unittest discover -s tests -t .. -v
```

These tests run entirely offline, using a mock hobrac/Mash setup and a loopback FASTA server, so they do not submit real Slurm jobs or download files from NCBI. They verify step overlap, temporary file cleanup, restarts, failure handling, and migration, and they include a 1,000-genome DAG dry run.
