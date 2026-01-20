import os

rule aln:
    input:
        reference = "reference/{accession}.fna",
        assembly = config["assembly"]
    output: "aln/vs_{accession}/aln.paf"
    threads: 12
    container: "docker://ghcr.io/cea-lbgb/hobrac-tools:latest"
    resources:
        mem_mb = config["minimap2_memory"],
        runtime = config["minimap2_runtime"]
    benchmark: "benchmarks/aln_{accession}.txt"
    shell: """
        minimap2 -x asm20 -t {threads} {input.reference} {input.assembly} > {output}
    """


rule gen_dgenies_index:
    input:
        aln = "aln/vs_{accession}/aln.paf",
        reference = "reference/{accession}.fna",
        assembly = config["assembly"]
    output: touch("aln/vs_{accession}/dgenies.done")
    params: 
        name = config["scientific_name"],
        assembly_prefix = config["scientific_name"].replace(" ", "_"),
        assembly_path = lambda wildcards, input: input.assembly if os.path.isabs(input.assembly) else f"../../{input.assembly}",
        reference_path = lambda wildcards, input: input.reference if os.path.isabs(input.reference) else f"../../{input.reference}"
    log: "logs/aln/gen_dgenies_index_{accession}.log"
    container: "docker://ghcr.io/cea-lbgb/hobrac-tools:latest"
    resources:
        mem_mb = 50000,
        runtime = 60
    benchmark: "benchmarks/dgenies_{accession}.txt"
    shell: """
        cd aln/vs_{wildcards.accession}

        dgenies_fasta_to_index -i {params.assembly_path} -n "{params.name}" -o query_{params.assembly_prefix}.idx
        dgenies_fasta_to_index -i {params.reference_path} -n "{wildcards.accession}" -o target_{wildcards.accession}.idx

        dotplotrs -m 2000 -p aln.paf -o dotplot_{params.assembly_prefix}_vs_{wildcards.accession}_significance.png --line-thickness 8
        dotplotrs -m 2000 -p aln.paf -o dotplot_{params.assembly_prefix}_vs_{wildcards.accession}_gravity.png --line-thickness 8 --gravity-ordering-only
        dotplotrs -m 2000 -p aln.paf -o dotplot_{params.assembly_prefix}_vs_{wildcards.accession}_bw.png --line-thickness 8 --no-color
        dotplotrs -m 2000 -p aln.paf -o dotplot_{params.assembly_prefix}_vs_{wildcards.accession}_bw_gravity.png --line-thickness 8 --no-color --gravity-ordering-only
    """


def get_all_ranking_targets(wildcards):
    checkpoint_output = checkpoints.select_references.get(**wildcards).output[0]
    accessions = []
    with open(checkpoint_output) as f:
        for line in f:
            if line.strip():
                accessions.append(line.strip())
    return expand(
        ["aln/vs_{accession}/dgenies.done", "aln/busco_{accession}"],
        accession=accessions
    )

rule rank_symlinks:
    input:
        selected_accessions = "mash/selected_accessions.txt",
        targets = get_all_ranking_targets
    output:
        "aln/ranking_symlinks.done"
    shell: """
        # Clean existing symlinks to avoid duplicates or stale links
        find aln/ -maxdepth 1 -type l -name "rank*_busco" -delete
        find aln/ -maxdepth 1 -type l -name "rank*_geno" -delete

        i=1
        while read -r accession; do
            # Create symlink for BUSCO
            if [ -d "aln/busco_${accession}" ]; then
                ln -s "busco_${accession}" "aln/rank${i}_busco"
            fi
            
            # Create symlink for genome-to-genome alignments
            if [ -d "aln/vs_${accession}" ]; then
                ln -s "vs_${accession}" "aln/rank${i}_geno"
            fi
            
            i=$((i+1))
        done < {input.selected_accessions}
        
        touch {output}
    """
