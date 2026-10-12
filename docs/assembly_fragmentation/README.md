# Assembly-fragmentation experiment

This workflow measures how hobrac behaves on fragmented assemblies. For each sampled ATLASea genome, it runs the complete hobrac pipeline three times: on the contigs, on the scaffolds and on the chromosomes of the same assembly. It then reports runtimes against N50, cumulative length and number of sequences, along with BUSCO dotplots arranged one genome per row.

```text
                                       ┌─ contigs ───── hobrac ─┐
Frozen manifest ── download ── levels ─┼─ scaffolds ─── hobrac ─┼─ report
                                       └─ chromosomes ─ hobrac ─┘
```

## Results and interpretation

The archived run evaluated **17 genomes from 17 phyla, 51 hobrac runs**, drawn with a fixed seed from the 386 ATLASea assemblies listed on 2 October 2026 and compared against MASH release 20260629. All 51 runs succeeded at the first attempt. The sample spans 7–3,505 sequences, 17 Mb–5.7 Gb and contig N50 from 1.2 to 30 Mb; six genomes have 222–1,194 scaffolds.

![Share of BUSCO anchors in significant pairs by assembly level, N50 and anchors per sequence](report/signal_vs_contiguity.png)

[Download the PDF](report/signal_vs_contiguity.pdf).

Each point is one run. The y-axis is the fraction of BUSCO anchors shared with the reference that fall in an assembly–reference sequence pair significant after Bonferroni correction, i.e. the anchors hobrac colours and uses. Colour is the phylum, shape the assembly level, and lines join the three levels of one genome.

### Synteny signal

**Scaffolds give the same signal as chromosomes in every genome.** The largest difference is 0.06 in *Acanthonyx lunulatus*, the one genome whose chromosome level selected a different reference; all other differences are at most 0.01. This includes the five genomes with hundreds of scaffolds: the unplaced scaffolds are short, carry few anchors, and neither add nor remove signal.

**Contigs lose signal only when they carry few anchors.** Thirteen genomes change by less than 0.10 between contigs and chromosomes. Four lose more, and they are the genomes whose contigs carry fewer than ten anchors on average:

| Genome | Contigs | Contig N50 | Anchors per anchor-bearing contig | Signal, contigs → scaffolds → chromosomes |
| --- | ---: | ---: | ---: | --- |
| *Ophelia neglecta* (Annelida) | 996 | 3.4 Mb | 5.4 | 0.08 → 0.51 → 0.51 |
| *Acanthonyx lunulatus* (Arthropoda) | 3,505 | 1.3 Mb | 2.6 | 0.32 → 0.82 → 0.75 |
| *Opheodesoma grisea* (Echinodermata) | 131 | 29.7 Mb | 9.4 | 0.45 → 0.72 → 0.72 |
| *Glandiceps talaboti* (Hemichordata) | 198 | 6.1 Mb | 5.9 | 0.77 → 0.93 → 0.93 |

No genome whose contigs carry more than 20 anchors on average lost signal. This is the expected behaviour of the test: a sequence pair is coloured only if the sequence holds enough anchors from a single reference chromosome to stay significant after correction over every assembly × reference pair, so the relevant quantity is the number of shared BUSCO genes per sequence, i.e. N50 relative to the BUSCO marker density of the lineage, rather than N50 alone. Short-range collinearity is still visible in grey in the contig dotplots of these genomes; what is lost is the chromosome-scale assignment.

**Four genomes have no significant pair at any level, chromosomes included.** Their limitation is the reference available in the MASH database, not fragmentation:

| Genome | Reference | MASH distance | Anchors | Likely cause |
| --- | --- | ---: | ---: | --- |
| *Rebecca billardiae* (Haptophyta) | *Prymnesium parvum* GCA_041296205.1 | 0.233 | 97 | Too few anchors (`eukaryota` lineage, 129 markers) |
| *Pseudobiceros splendidus* (Platyhelminthes) | *Schmidtea mediterranea* GCA_045838265.1 | 0.238 | 944 | Polyclad vs triclad; extensively rearranged reference |
| *Dudresnaya verticillata* (Rhodophyta) | *Pyropia yezoensis* GCA_009829735.1 | 0.291 | 1,072 | Florideophyceae vs Bangiophyceae |
| *Halophila stipulacea* (Streptophyta) | *Allium cepa* GCA_030765085.1 | 0.195 | 2,251 | Alismatales vs Asparagales; polyploid history, 16 Gb reference |

MASH distance alone does not predict this: *Ophelia neglecta* and *Pseudobiceros splendidus* both sit at 0.238 and only the former works. These four rows are kept in all figures and tables; they show the method's dependence on the reference, which the [MASH reference survey](../mash_reference_survey/README.md) measures across phyla.

### Reference selection

The three levels of a genome selected the same reference in 16 of 17 genomes (`reference_selection.tsv`). In *Acanthonyx lunulatus* the chromosome level selected GCA_054960195.1 instead of GCA_023348025.1, at a MASH distance of 0.1613 versus 0.1610. Removing the unplaced scaffolds changed the sketch enough to swap two references that are equally close; the signal was similar with both.

### Runtime and memory

![hobrac runtime per step against the number of sequences](report/runtime_vs_sequences.png)

[Download the PDF](report/runtime_vs_sequences.pdf); the same figure against [N50](report/runtime_vs_n50.pdf) and [cumulative length](report/runtime_vs_total_length.pdf).

End-to-end wall time ranged from 3.1 min (*Dudresnaya verticillata*, 76 Mb) to 161 min (*Halophila stipulacea*, 5.7 Gb) on 24 CPUs of one node type. **The number of sequences has no measurable effect.** The contig level ran in 0.81–1.09× the time of the chromosome level, the scaffold level in 0.83–1.28×; the two largest deviations are the reference change in *Acanthonyx lunulatus* and one slower scaffold run of *Striarca lactea*. The longest BUSCO run (assembly or reference, which execute concurrently) accounts for 51–95% of wall time and depends on genome size and lineage, not on contiguity. Minimap2 never exceeded 28% of a run and was *faster* on contigs than on chromosomes for the largest genomes (*Halophila stipulacea* 12 vs 44 min, *Striarca lactea* 4.5 vs 18 min, *Pseudobiceros splendidus* 1.0 vs 3.5 min). Peak memory (`benchmarks.tsv`, max RSS over all rules) was 2.1–130 GB; it is set by BUSCO on the larger of the two genomes, the extreme being the 16 Gb *Allium cepa* reference, and the contig level used 0.65–1.20× the memory of the chromosome level.

### Conclusion

On these 17 genomes, **hobrac is insensitive to the number of sequences and loses signal only when individual sequences carry too few shared BUSCO genes**, which happened at the contig level for N50 below a few megabases in gene-sparse genomes. Scaffold-level assemblies with hundreds of unplaced scaffolds behave like their chromosomes, in signal and in runtime.

The experiment uses the three levels that exist for every real chromosome-level assembly rather than artificial breakpoints, so the fragmentation of each genome is what its assembler produced, not a controlled gradient. It gives three points per genome; it does not locate the exact N50 at which signal disappears, and the dotplot significance is the only output quantified. The four genomes without signal are limited by their reference and say nothing about fragmentation.

### Archived inputs and outputs

`report/` holds the output of `summarize` for the archived run: `runs.tsv`, `benchmarks.tsv`, `reference_selection.tsv`, `summary.json`, the figures and the dotplot grids (see [Outputs](#outputs)). It is a **report archive**, not a runnable work directory: the frozen manifest, taxdump, downloaded genomes, level FASTAs and hobrac run directories are not included. The runtime and signal figures were regenerated from `runs.tsv` after the run with a revised layout; `script_sha256` in `summary.json` identifies the script of the original run, not the later plotting changes.

## Design

### Genome set

Preparation downloads the NCBI Datasets reports for every assembly in the ATLASea umbrella BioProject (`PRJEB64126`). It writes them to `atlasea_assemblies.tsv` with their phylum and an eligibility status. Phyla come from a frozen NCBI taxdump, which hobrac's taxonkit then reuses during the runs. An assembly is excluded when:

- it is a RefSeq copy (`GCF_`), an alternate haplotype, not current, or below chromosome level,
- its taxid is missing from the frozen taxdump, or its lineage has no phylum or species rank,
- its phylum has no sketch in the frozen MASH release, so hobrac could not select a reference automatically,
- another assembly represents the same species. The representative is chosen by assembly level (Complete Genome first), then reference-genome status, then newest release, then accession.

One species per phylum is then drawn with a fixed seed. Pools are sorted by accession, so the draw does not depend on the order of the API results. On 2 October 2026, the BioProject contained 386 assemblies, and the draw selected 17 genomes from 17 phyla, which means 51 hobrac runs.

### Assembly levels

Each level is built from the NCBI FASTA and its sequence report.

| Level | Sequences | Names |
| --- | --- | --- |
| `chromosomes` | Assembled molecules | `chr<name>` |
| `scaffolds` | Chromosomes, plus unlocalized and unplaced scaffolds | GenBank accessions for non-chromosomes |
| `contigs` | Every scaffold-level sequence split at runs of at least 10 Ns | `<name>_ctg<i>` when split |

`levels.json` records the number of sequences, cumulative length, N50 and remaining Ns of each level. It also stores NCBI's contig count, which can be compared with the `contigs` level.

### hobrac runs

Each genome × level runs the default pipeline end to end: automatic reference selection with MASH, BUSCO on the assembly and on the reference, the minimap2 alignment, and the plots:

```text
hobrac -a <level>.fna.gz -n <name> -t <taxid> -o runs/<accession>/<level>/hobrac \
    --mash-db-version <frozen release> -e local [--use-<container>] [HOBRAC_ARGS]
```

- **Reference selection.** The reference is not fixed. The MASH release (`latest` by default) is resolved once during preparation, and every run uses that dated release, so a difference in the selected reference comes from the assembly level alone. `reference_selection.tsv` reports whether the three levels of each genome selected the same reference.
- **One shared allocation for the entire run.** Each run shares a Slurm allocation (`CPUS`, `MEMORY`, `TIME`), and hobrac executes all of its steps inside it. The end-to-end time therefore excludes queueing. hobrac uses every CPU that Slurm assigns to the job.
- **Fixed options.** `HOBRAC_ARGS` applies the same extra options to every run. Options that would change the experiment are rejected: input, output, reference, executor, MASH release, precomputed BUSCO, `--stop-after-mash`, `--skip-genomic` and `--ref-count`.
- **Retries.** A failed run is restarted from scratch, so every reported runtime comes from a single uninterrupted execution. `attempt` in `runs.tsv` counts the tries.
- **Completed runs.** A successful run is never recomputed. Delete `runs/<accession>/<level>/` to force a rerun.

## Outputs

All reports are written to `$WORK/report/`:

| File | Content |
| --- | --- |
| `runs.tsv` | One row per genome × level. It records the status, sequences, cumulative length and N50, as well as the selected reference and its MASH distance. It also includes the BUSCO lineage and C/S/D/F/M scores, BUSCO anchors and the sequences carrying them, and the share of anchors in significant dotplot pairs. Finally, it lists the end-to-end and per-step-group wall times, the node, the CPU model and the number of CPUs. |
| `benchmarks.tsv` | Every hobrac Snakemake benchmark, with wall time, CPU time and maximum RSS. |
| `reference_selection.tsv` | Selected reference and distance for each level, and `same_reference` (`yes`, `no` or `incomplete`). |
| `runtime_vs_n50.*`, `runtime_vs_total_length.*`, `runtime_vs_sequences.*` | Runtime per step against each variable on log–log axes. Lines join the levels of one genome. |
| `signal_vs_contiguity.*` | Share of BUSCO anchors in significant dotplot pairs, against assembly level, N50 and anchors per anchor-bearing sequence. Lines join the levels of one genome. |
| `dotplot_grid_busco.pdf`, `dotplot_grid_busco_bw.pdf` | hobrac's BUSCO dotplots (`aln/busco_<reference>/dotplot_busco*.png`), with one genome per row and the contig, scaffold and chromosome levels as columns. Pages are also saved as PNGs in `dotplot_grids/`. |
| `summary.json` | Status counts, reference agreement and provenance hashes. |

Each runtime figure has one panel per step group, built from hobrac's benchmark files:

| Panel | hobrac rules |
| --- | --- |
| End-to-end hobrac run | the whole `hobrac` command |
| BUSCO (assembly) | `busco_assembly` |
| BUSCO (reference) | `busco_reference_*` |
| Minimap2 alignment | `aln_*` |
| MASH search | `mash`, `select_references` |
| Downloads | `download_db`, `get_reference_*`, `get_busco_dataset`, `download_busco_dataset` |
| Other steps | all remaining rules |

Notes on interpretation:

- **Downloads** depend on the network, not on the assembly.
- **BUSCO on the reference** reprocesses identical input only when hobrac selects the same reference at all three levels of a genome (`same_reference` in `reference_selection.tsv`). For those genomes, the spread of these times shows the timing noise of the cluster. When the selected reference changes, the times also reflect a different reference genome.
- **Node type.** Runtimes are comparable only on one node type. Constrain the jobs with `EXTRA_SBATCH`; `runs.tsv` records the node and CPU model of every run.
- **Colours in the BUSCO dotplots.** dotplotrs colours an assembly–reference sequence pair only when its Fisher test, Bonferroni-corrected over the number of assembly × reference sequences, is significant. Other anchors are drawn in grey. Sequences without a complete BUSCO shared with the reference are not drawn at all.

## Running the workflow

You need Linux, Bash, Python >=3.11, installed hobrac with its tools (or a container runtime), Snakemake 9 with its Slurm executor plugin, and matplotlib. Preparation and the hobrac runs need network access to NCBI and to the Genoscope MASH server. If needed, you can install the workflow dependencies directly into your hobrac environment:

```bash
python3 -m pip install 'snakemake>=9,<10' snakemake-executor-plugin-slurm 'matplotlib>=3.6'
```

Copy **this whole directory** to your machine, and make sure the copy includes `scripts/` and `profiles/`.

### 1. Configure

```bash
cp config.example.sh config.local.sh
# Edit config.local.sh before submitting anything.
```

Set `WORK` to an absolute shared path without spaces or shell metacharacters. Then fill in your environment, the Slurm settings (including a node constraint in `EXTRA_SBATCH`) and the per-run resources `CPUS`, `MEMORY` and `TIME`. Both BUSCO steps request 12 threads, so with 24 CPUs they can run at the same time. `MAX_CONCURRENT` limits the submitted child jobs, and `MASTER_TIME` must cover the slowest hobrac run plus queueing.

**A different seed, BioProject or MASH release requires a new `WORK` directory.**

### 2. Prepare and review

```bash
bash submit.sh prepare
```

Before you proceed, inspect `manifest.tsv`, `phylum_coverage.tsv`, `atlasea_assemblies.tsv` and `sampling_exclusions.json`.

### 3. Run

```bash
bash submit.sh workflow
```

By default (`CLEANUP=true`), the downloaded and level FASTAs are deleted once their hobrac runs succeed. hobrac already deletes its own FASTA copies and BUSCO lineage downloads.

### 4. Inspect results

Start with `$WORK/report/runs.tsv`, the runtime figures and the dotplot grids.

## Restart or inspect a run

```bash
bash submit.sh workflow                  # Continue unfinished work
bash submit.sh summarize                 # Produce a partial report separately
bash workflow.sh config.local.sh --dry-run
```

To run the workflow locally instead of through Slurm, for example inside one allocation so that every run executes serially on the same node:

```bash
bash workflow.sh config.local.sh --executor local
```

The local executor uses `profiles/local` and `--cores CPUS`, so exactly one hobrac run executes at a time (`rule hobrac` declares `threads: CPUS`); `MAX_CONCURRENT` does not apply. Append `--cores 2×CPUS` to run two at once on a larger machine, but note that each hobrac run then still sees every CPU of the allocation. Never pass `--jobs` to the local executor: Snakemake stops restricting cores as soon as a node limit is set, and `MAX_CONCURRENT` runs would start together.

A serial series on one dedicated node, with the same per-run CPU count as the Slurm runs:

```bash
#!/bin/bash
#SBATCH --nodes=1 --exclusive --cpus-per-task=24 --mem=160G --time=2-00:00:00
bash /path/to/assembly_fragmentation/workflow.sh config.local.sh --executor local
```

`--cpus-per-task` must equal `CPUS`; the inner hobrac uses every CPU in its affinity mask, and `runs.tsv` records that count as `cpus`. Expect roughly one day for the 51 runs of the default sample on a current Xeon node.

Never run two controllers against the same `WORK` directory.

## Directory layout

| Path | Purpose |
| --- | --- |
| `config.example.sh` | Template configuration to copy and edit for your run |
| `submit.sh` | Submits preparation, the controller, or reporting jobs to Slurm |
| `workflow.sh`, `Snakefile` | Launches the Snakemake controller and defines the DAG |
| `profiles/` | Scheduler profiles: `slurm` for the default controller, `local` for `--executor local` |
| `scripts/` | Python backend and Slurm job entry point |
| `tests/` | Offline unit tests |
| `report/` | Archived result tables, figures and dotplot grids of the run described above |

To run the test suite from this directory:

```bash
python3 -m unittest discover -s tests -t .. -v
```
