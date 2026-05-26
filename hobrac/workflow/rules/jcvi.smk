def get_busco_reference_dirs(wildcards):
    """Get all BUSCO reference directories based on selected accessions."""
    checkpoint_output = checkpoints.select_references.get(**wildcards).output[0]
    accessions = []
    with open(checkpoint_output) as f:
        for line in f:
            if line.strip():
                accessions.append(line.strip())
    return expand("busco/busco_reference_{accession}", accession=accessions)


rule resolve_jcvi_color_scheme:
    input:
        chosen_dataset="busco/chosen_dataset.txt",
    output:
        resolved="aln/jcvi_karyotype/resolved_colors.txt",
    params:
        scheme=config.get("jcvi_color_scheme", ""),
    run:
        import os, sys

        resolved_path = ""
        scheme = params.scheme
        if scheme:
            with open(input.chosen_dataset) as f:
                dataset = f.readline().strip().split("\t")[1]
            colors_dir = os.path.join(
                os.path.dirname(os.path.dirname(workflow.snakefile)),
                "colors",
            )
            color_file = os.path.join(
                colors_dir, f"Busco.Colors.{scheme}.{dataset}"
            )
            if os.path.isfile(color_file):
                resolved_path = color_file
            else:
                logger.warning(
                    f"No pre-computed {scheme} color file for dataset"
                    f" '{dataset}'. Falling back to default coloring."
                )
        with open(output.resolved, "w") as out:
            out.write(resolved_path)


rule jcvi_synteny:
    input:
        busco_assembly="busco/busco_assembly",
        busco_references=get_busco_reference_dirs,
        accession_order="mash/selected_accessions.txt",
        assembly=config["assembly"],
        resolved_colors="aln/jcvi_karyotype/resolved_colors.txt",
    output:
        seqids="aln/jcvi_karyotype/seqids",
        layouts="aln/jcvi_karyotype/layouts",
    benchmark:
        "benchmarks/jcvi_synteny.txt"
    resources:
        mem_mb=8000,
        runtime=30,
    params:
        manual_refs=config.get("manual_references", ""),
        assembly_name=config["scientific_name"].replace(" ", "_"),
        outdir="aln/jcvi_karyotype",
        min_busco_genes=config.get("min_busco_genes", 0),
        jcvi_custom_colors=config.get("jcvi_custom_colors", ""),
        jcvi_names=config.get("jcvi_names", ""),
        hide_non_significant=config.get("hide_non_significant", False),
        skip_alg=config.get("skip_alg", False),
        jcvi_pvalue=config.get("jcvi_pvalue", 0.01),
        jcvi_min_chain_genes=config.get("jcvi_min_chain_genes", 5),
    shell:
        """
        RESOLVED_COLORS=$(cat {input.resolved_colors})
        if [ -n "{params.jcvi_custom_colors}" ]; then
            COLOR_ARG="{params.jcvi_custom_colors}"
        elif [ -n "$RESOLVED_COLORS" ]; then
            COLOR_ARG="$RESOLVED_COLORS"
        else
            COLOR_ARG=""
        fi

        jcvi_synteny \
            --busco_assembly {input.busco_assembly}/run*/full_table.tsv \
            --assembly-fasta {input.assembly} \
            --busco_references {input.busco_references} \
            --accession_order {input.accession_order} \
            --manual_refs "{params.manual_refs}" \
            --assembly_name "{params.assembly_name}" \
            --output_dir {params.outdir} \
            --min-busco-genes {params.min_busco_genes} \
            --jcvi-pvalue {params.jcvi_pvalue} \
            --jcvi-min-chain-genes {params.jcvi_min_chain_genes} \
            $([ -n "$COLOR_ARG" ] && echo "--jcvi-custom-colors $COLOR_ARG") \
            $([ -n "{params.jcvi_names}" ] && echo "--jcvi-names {params.jcvi_names}") \
            $([ "{params.hide_non_significant}" = "True" ] && echo "--hide-non-significant") \
            $([ "{params.skip_alg}" = "True" ] && echo "--skip-alg")
        """


rule jcvi_karyotype:
    input:
        seqids="aln/jcvi_karyotype/seqids",
        layouts="aln/jcvi_karyotype/layouts",
    output:
        "aln/jcvi_karyotype/karyotype.png",
    benchmark:
        "benchmarks/jcvi_karyotype.txt"
    container:
        HOBRAC_TOOLS
    resources:
        mem_mb=4000,
        runtime=10,
    shell:
        """
        cd aln/jcvi_karyotype
        python -m jcvi.graphics.karyotype seqids layouts \
            --dpi 100 --figsize 12x10 --notex --basepair -o karyotype.png
        """
