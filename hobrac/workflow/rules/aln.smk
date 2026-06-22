import os


rule aln:
    input:
        reference="reference/{accession}.fna",
        assembly=config["assembly"],
    output:
        "aln/vs_{accession}/aln.paf",
    benchmark:
        "benchmarks/aln_{accession}.txt"
    container:
        HOBRAC_TOOLS
    threads: 12
    resources:
        mem_mb=config["minimap2_memory"],
        runtime=config["minimap2_runtime"],
    shell:
        """
        minimap2 -x asm20 -t {threads} {input.reference} {input.assembly} > {output}
    """


rule gen_dgenies_index:
    input:
        aln="aln/vs_{accession}/aln.paf",
        reference="reference/{accession}.fna",
        assembly=config["assembly"],
    output:
        touch("aln/vs_{accession}/dgenies.done"),
    log:
        "logs/aln/gen_dgenies_index_{accession}.log",
    benchmark:
        "benchmarks/dgenies_{accession}.txt"
    container:
        HOBRAC_TOOLS
    resources:
        mem_mb=50000,
        runtime=60,
    params:
        name=config["scientific_name"],
        assembly_prefix=config["scientific_name"].replace(" ", "_"),
        assembly_path=lambda wildcards, input: (
            input.assembly
            if os.path.isabs(input.assembly)
            else f"../../{input.assembly}"
        ),
        reference_path=lambda wildcards, input: (
            input.reference
            if os.path.isabs(input.reference)
            else f"../../{input.reference}"
        ),
    shell:
        """
        cd aln/vs_{wildcards.accession}

        dgenies_fasta_to_index -i {params.assembly_path} -n "{params.name}" -o query_{params.assembly_prefix}.idx
        dgenies_fasta_to_index -i {params.reference_path} -n "{wildcards.accession}" -o target_{wildcards.accession}.idx

        dotplotrs -m 2000 -p aln.paf -o dotplot.png --line-thickness 4 --dump-significance significance_results.txt
        dotplotrs -m 2000 -p aln.paf -o dotplot_bw.png --line-thickness 4 --no-color
    """


def get_all_ranking_targets(wildcards):
    checkpoint_output = checkpoints.select_references.get(**wildcards).output[0]
    accessions = []
    with open(checkpoint_output) as f:
        for line in f:
            if line.strip():
                accessions.append(line.strip())
    patterns = ["aln/busco_{accession}"]
    if not config.get("skip_genomic", False):
        patterns.append("aln/vs_{accession}/dgenies.done")
    return expand(patterns, accession=accessions)


def get_ranking_report_inputs(wildcards):
    """NCBI assembly reports, so the ranked symlinks can carry organism names.

    Each report is co-produced with its reference .fna by get_reference, so it
    already exists by ranking time; declaring it makes that dependency explicit.
    Manual references have no NCBI report, so request none (requiring one would
    force a download) -- rank_symlinks then falls back to the accession.
    """
    if config.get("manual_references"):
        return []
    checkpoint_output = checkpoints.select_references.get(**wildcards).output[0]
    accessions = []
    with open(checkpoint_output) as f:
        for line in f:
            if line.strip():
                accessions.append(line.strip())
    return expand(
        "reference/{accession}_assembly_report.txt", accession=accessions
    )


rule rank_symlinks:
    input:
        selected_accessions="mash/selected_accessions.txt",
        targets=get_all_ranking_targets,
        reports=get_ranking_report_inputs,
    output:
        "aln/ranking_symlinks.done",
    run:
        import glob
        import os

        from hobrac.jcvi_synteny.io import parse_organism_name

        # Clean existing symlinks (old or new naming) to avoid stale links.
        # The rank* glob matches both "rank1_busco" and "rank1_<Species>_busco".
        for pattern in ("aln/rank*_busco", "aln/rank*_geno"):
            for link in glob.glob(pattern):
                if os.path.islink(link):
                    os.unlink(link)

        with open(input.selected_accessions) as f:
            accessions = [line.strip() for line in f if line.strip()]

        # The rank<i> prefix keeps every symlink unique even when two genomes
        # share a species, so identical organism labels never collide.
        for i, accession in enumerate(accessions, start=1):
            organism = parse_organism_name(
                f"reference/{accession}_assembly_report.txt"
            )
            label = organism.replace(" ", "_") if organism else accession
            for src, suffix in (
                (f"busco_{accession}", "busco"),
                (f"vs_{accession}", "geno"),
            ):
                if os.path.isdir(os.path.join("aln", src)):
                    os.symlink(src, f"aln/rank{i}_{label}_{suffix}")

        open(output[0], "w").close()
