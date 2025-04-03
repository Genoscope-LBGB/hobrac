rule aln:
    input:
        mash_output = rules.select_closest_reference.output,
        assembly = config["assembly"]
    output: touch("aln/aln.done")
    threads: 12
    resources:
        mem_mb = 60000,
        runtime = 4 * 60
    benchmark: "benchmarks/aln.txt"
    shell: """
        prefix=$(cat {input.mash_output} | cut -f 1)
        filepath=$(cat {input.mash_output} | cut -f 2)

        mkdir -p aln/vs_${{prefix}} && cd aln/vs_${{prefix}}

        minimap2 -x asm20 -t {threads} ${{filepath}} {input.assembly} > aln.paf
    """


rule gen_dgenies_index:
    input:
        rules.aln.output,
        mash_output = rules.select_closest_reference.output,
        assembly = config["assembly"]
    output: touch("aln/get_dgenies_index.done")
    params: 
        name = config["scientific_name"],
        assembly_prefix = config["scientific_name"].replace(" ", "_")
    resources:
        mem_mb = 10000,
        runtime = 60
    benchmark: "benchmarks/aln.txt"
    shell: """
        ref_prefix=$(cat {input.mash_output} | cut -f 1)
        ref_name=$(echo $ref_prefix | sed 's/_/ /g')
        ref_filepath=$(cat {input.mash_output} | cut -f 2)

        cd aln/vs_${{ref_prefix}}

        dgenies_fasta_to_index -i {input.assembly} -n "{params.name}" -o query_{params.assembly_prefix}.idx
        dgenies_fasta_to_index -i ${{ref_filepath}} -n "${{ref_name}}" -o target_${{ref_prefix}}.idx

        
        dotplot -p aln.paf -o {output}/dotplot_{params.assembly_prefix}_vs_${{ref_prefix}}.png
    """
