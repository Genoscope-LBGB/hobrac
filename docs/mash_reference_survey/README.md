# Mash reference-availability survey

## What this is

We want to know how far away the nearest eligible reference tends to be when hobrac stops after the Mash step, across a broad set of eukaryotic species. The approach: sample 100 species balanced across the phyla that hobrac supports, download one NCBI assembly per species, run the installed hobrac with `--stop-after-mash` only, and plot the nearest eligible reference distance grouped by query phylum. The default seed is `20260717`.

## 1. Configure your Slurm cluster run

Copy this whole directory to the cluster:

```bash
cd docs/mash_reference_survey
cp config.example.sh config.local.sh
# Edit config.local.sh before submitting anything.
```

The settings that matter most:

**`WORK`** should be an absolute path, with no spaces or shell metacharacters. You need room for 100 compressed genomes plus hobrac's uncompressed prepared FASTAs, the taxonomy, sketches and logs. Large eukaryotic genomes can push the total to hundreds of GB or more, and nothing is deleted automatically.

**`ENV_SETUP`** is optional but useful: the absolute path of a Bash script that activates your existing hobrac environment. It runs both at submission time and inside each job.

**`PARTITION`, `ACCOUNT`, `QOS`** take your site's values, or leave them empty for defaults.

**`CONTAINER`** should be `none` if your dependencies are native, or `apptainer`/`singularity` if that is how your installed hobrac normally runs.

**Resources** default to 2 CPUs, 24 GB and 4 hours per query, with at most 4 concurrent queries. Preparation and reporting each request 16 GB because the frozen NCBI taxonomy is loaded in memory. Adjust for your cluster and the genomes you expect.

### Scope decision

`ANCESTOR_TAXID=2759` samples across supported eukaryotic phyla, including animals, fungi, plants and protists. Use `33208` instead if the intended article comparison should be animals only. Decide this before you look at any distances.

### Required software and access

You will need Linux, Bash, Slurm (`sbatch`).

Python >=3.11, your installed `hobrac` with `--stop-after-mash` and `-e local` support, and host-side `mash` (including for container runs, since preparation validates downloaded sketches with `mash info` and execution records the host version).

Matplotlib for the report: install it in your chosen environment if it is not already there, e.g. `python3 -m pip install 'matplotlib>=3.6'`.

`job.sh` sets `TAXONKIT_DB` to the frozen survey taxonomy for query jobs. For Apptainer/Singularity, the runner appends an explicit bind of the canonical `WORK` path to the container's bind environment, preserving whatever binds your site already configures. This makes the shared sketch symlink target visible under the same absolute path inside the container. Confirm your site's container policy permits this during the pilot. No nested Slurm submission occurs: each task calls `hobrac ... -e local --stop-after-mash` directly.

## 2. Freeze metadata and select the queries

```bash
bash submit.sh prepare
```

Wait for the preparation job to finish successfully. Its Slurm log is under `$WORK/logs/`. Then inspect:

- `manifest.tsv`: the exact query accessions with versions, taxids, species taxids, names, phyla, assembly levels, release dates, download URLs and zero-based task IDs.
- `phylum_coverage.tsv`: eligible assembly and species counts, sketch availability, and allocated sample sizes. Phyla with no eligible queries are retained here so you can see the gaps.
- `sampling_exclusions.json`: counts excluded by each sequential metadata filter.

Preparation freezes NCBI metadata, taxonomy, and one hobrac sketch per sampled phylum, rather than downloading a potentially different reference database for every query. The `snapshot.json` hashes these inputs. **Never regenerate a finished snapshot in place.** A different scope, seed, or fresh database requires a new `WORK` directory. Re-running an interrupted preparation reuses already-downloaded files after checksum verification; a completed preparation refuses to overwrite its snapshot.

### How the sampling works

1. Download GenBank assembly summaries for `invertebrate`, `vertebrate_mammalian`, `vertebrate_other`, `plant`, `fungi` and `protozoa`.

2. Keep latest Chromosome / Complete Genome, full-genome assemblies with no `excluded_from_refseq` reason.

3. Require a species ancestor and a named phylum in the frozen NCBI taxonomy, species-taxid agreement with the assembly summary, ancestry within the selected scope, a downloadable versioned GCA accession, and an exact matching phylum sketch in hobrac's published directory.

4. Choose one assembly per species deterministically: Complete Genome before Chromosome, then reference/representative designation, then newest release date, then lexicographic accession as the final tie-break.

5. Seeded random sampling without replacement within each phylum. Allocate one query to each phylum per round in a seeded random phylum order until the total is reached. Phyla with too few eligible species contribute all available species; their unused allocation is redistributed. If fewer species than the requested total exist, preparation fails rather than silently reducing the count.

`eligible_species.tsv` preserves the entire sampling frame.

## 3. Run the Slurm array

```bash
bash submit.sh array
```

This submits one task per manifest row and an `afterany` summary job and the report runs even if some tasks fail. 

Each task downloads its query with bounded retries, validates gzip integrity, and quarantines invalid gzip or FASTA responses under `rejected_downloads/` before retrying the same accession (at most three validation attempts). It then links the frozen phylum sketch into hobrac's expected `mash/mash_db.msh` path, and stops after producing `mash/mash.dist`. All raw reference comparisons are kept. Current hobrac uses `mash dist -s 10000` and the reference sketch's stored k-mer parameters, so inspect the saved `sketches/*.info.txt` for the actual sketch settings.

## 4. Inspect the results

Under `$WORK/report/`:

| File | Content |
| --- | --- |
| `nearest_references.tsv` | One row for every sampled query, including failures: nearest eligible reference, distance, P-value, matching hashes, taxonomic relationship and exclusion counts |
| `distance_by_phylum.png`, `.pdf` | Individual distances and boxplots by query phylum, with sample counts |
| `phylum_summary.tsv` | Counts of sampled queries, available references, no-shared-hash cases and incomplete/failed queries |
| `summary.json` | Overall status counts and provenance |
| `retry_task_ids.txt` | Comma-separated ordinary failed/incomplete task IDs for retry |
| `reset_required_task_ids.txt` | Corrupted/invalid results needing investigation before rerun |