rule aln:
    input:
        mash_output = rules.select_closest_reference.output,
        reference = rules.get_top_reference.output,
        assembly = config["assembly"]
    output: touch("aln/aln.done")
    log: "logs/aln/aln.log"
    threads: 12
    resources:
        mem_mb = config["minimap2_memory"],
        runtime = 12 * 60
    benchmark: "benchmarks/aln.txt"
    shell: """
        prefix=$(cat {input.mash_output})
        filepath=$(cat {input.reference})

        mkdir -p aln/vs_${{prefix}} && cd aln/vs_${{prefix}}

        minimap2 -x asm20 -t {threads} ${{filepath}} {input.assembly} > aln.paf
    """


rule gen_dgenies_index:
    input:
        rules.aln.output,
        mash_output = rules.select_closest_reference.output,
        reference = rules.get_top_reference.output,
        assembly = config["assembly"]
    output: touch("aln/get_dgenies_index.done")
    params: 
        name = config["scientific_name"],
        assembly_prefix = config["scientific_name"].replace(" ", "_")
    log: "logs/aln/gen_dgenies_index.log"
    resources:
        mem_mb = 10000,
        runtime = 60
    benchmark: "benchmarks/aln.txt"
    shell: """
        accession=$(cat {input.mash_output})
        ref_filepath=$(cat {input.reference})

        cd aln/vs_${{accession}}

        dgenies_fasta_to_index -i {input.assembly} -n "{params.name}" -o query_{params.assembly_prefix}.idx
        dgenies_fasta_to_index -i ${{ref_filepath}} -n "${{accession}}" -o target_${{accession}}.idx

        dotplotrs -m 2000 -p aln.paf -o dotplot_{params.assembly_prefix}_vs_${{accession}}.png --line-thickness 2
    """
