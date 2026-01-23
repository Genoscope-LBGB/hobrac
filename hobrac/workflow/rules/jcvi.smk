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
        accession_order = "mash/selected_accessions.txt"
    output:
        directory("aln/jcvi_karyotype")
    params:
        manual_refs = config.get("manual_references", ""),
        assembly_name = config["scientific_name"].replace(" ", "_")
    container: "docker://ghcr.io/cea-lbgb/hobrac-tools:latest"
    resources:
        mem_mb = 8000,
        runtime = 30
    benchmark: "benchmarks/jcvi_synteny.txt"
    shell:
        """
        jcvi_synteny \
            --busco_assembly {input.busco_assembly}/run*/full_table.tsv \
            --busco_references {input.busco_references} \
            --accession_order {input.accession_order} \
            --manual_refs "{params.manual_refs}" \
            --assembly_name "{params.assembly_name}" \
            --output_dir {output}
        """


rule jcvi_karyotype:
    input:
        "aln/jcvi_karyotype"
    output:
        "aln/jcvi_karyotype/karyotype.png"
    container: "docker://ghcr.io/cea-lbgb/hobrac-tools:latest"
    resources:
        mem_mb = 4000,
        runtime = 10
    benchmark: "benchmarks/jcvi_karyotype.txt"
    shell:
        """
        cd {input}
        python -m jcvi.graphics.karyotype seqids layouts \
            --dpi 100 --figsize 12x10 --basepair -o karyotype.png
        """
