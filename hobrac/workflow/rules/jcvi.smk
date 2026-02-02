def get_busco_reference_dirs(wildcards):
    """Get all BUSCO reference directories based on selected accessions."""
    checkpoint_output = checkpoints.select_references.get(**wildcards).output[0]
    accessions = []
    with open(checkpoint_output) as f:
        for line in f:
            if line.strip():
                accessions.append(line.strip())
    return expand("busco/busco_reference_{accession}", accession=accessions)


rule jcvi_synteny:
    input:
        busco_assembly = "busco/busco_assembly",
        busco_references = get_busco_reference_dirs,
        accession_order = "mash/selected_accessions.txt",
        assembly = config["assembly"]
    output:
        seqids = "aln/jcvi_karyotype/seqids",
        layouts = "aln/jcvi_karyotype/layouts"
    params:
        manual_refs = config.get("manual_references", ""),
        assembly_name = config["scientific_name"].replace(" ", "_"),
        outdir = "aln/jcvi_karyotype",
        min_busco_genes = config.get("min_busco_genes", 0),
        jcvi_custom_colors = config.get("jcvi_custom_colors", ""),
        jcvi_names = config.get("jcvi_names", ""),
        hide_non_significant = config.get("hide_non_significant", False)
    resources:
        mem_mb = 8000,
        runtime = 30
    benchmark: "benchmarks/jcvi_synteny.txt"
    shell:
        """
        jcvi_synteny \
            --busco_assembly {input.busco_assembly}/run*/full_table.tsv \
            --assembly-fasta {input.assembly} \
            --busco_references {input.busco_references} \
            --accession_order {input.accession_order} \
            --manual_refs "{params.manual_refs}" \
            --assembly_name "{params.assembly_name}" \
            --output_dir {params.outdir} \
            --min-busco-genes {params.min_busco_genes} \
            $([ -n "{params.jcvi_custom_colors}" ] && echo "--jcvi-custom-colors {params.jcvi_custom_colors}") \
            $([ -n "{params.jcvi_names}" ] && echo "--jcvi-names {params.jcvi_names}") \
            $([ "{params.hide_non_significant}" = "True" ] && echo "--hide-non-significant")
        """


rule jcvi_karyotype:
    input:
        seqids = "aln/jcvi_karyotype/seqids",
        layouts = "aln/jcvi_karyotype/layouts"
    output:
        "aln/jcvi_karyotype/karyotype.png"
    container: HOBRAC_TOOLS
    resources:
        mem_mb = 4000,
        runtime = 10
    benchmark: "benchmarks/jcvi_karyotype.txt"
    shell:
        """
        cd aln/jcvi_karyotype
        python -m jcvi.graphics.karyotype seqids layouts \
            --dpi 100 --figsize 12x10 --notex --basepair -o karyotype.png
        """
