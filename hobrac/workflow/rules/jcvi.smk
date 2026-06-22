def get_busco_reference_dirs(wildcards):
    """Get all BUSCO reference directories based on selected accessions."""
    checkpoint_output = checkpoints.select_references.get(**wildcards).output[0]
    accessions = []
    with open(checkpoint_output) as f:
        for line in f:
            if line.strip():
                accessions.append(line.strip())
    return expand("busco/busco_reference_{accession}", accession=accessions)


def get_dotplot_grid_inputs(wildcards):
    """Per-reference dotplots + assembly reports feeding the dotplot grid."""
    checkpoint_output = checkpoints.select_references.get(**wildcards).output[0]
    accessions = []
    with open(checkpoint_output) as f:
        for line in f:
            if line.strip():
                accessions.append(line.strip())
    patterns = [
        "aln/jcvi_karyotype/dotplots/{accession}.png",
        "aln/jcvi_karyotype/dotplots/{accession}_dark.png",
    ]
    # Manual references have no NCBI assembly report; only get_reference can
    # produce one (by downloading), so requiring it would re-trigger a download.
    # grid.py falls back to --jcvi-names / accession when the report is absent.
    if not config.get("manual_references"):
        patterns.append("reference/{accession}_assembly_report.txt")
    return expand(patterns, accession=accessions)


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
            import hobrac

            with open(input.chosen_dataset) as f:
                dataset = f.readline().strip().split("\t")[1]
            colors_dir = os.path.join(
                os.path.dirname(hobrac.__file__),
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
        gene_chains="aln/jcvi_karyotype/gene_chains.tsv",
        orders=directory("aln/jcvi_karyotype/dotplot_orders"),
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
        jcvi_permissive_alg=config.get("jcvi_permissive_alg", False),
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
            $([ "{params.skip_alg}" = "True" ] && echo "--skip-alg") \
            $([ "{params.jcvi_permissive_alg}" = "True" ] && echo "--jcvi-permissive-alg")
        """


rule jcvi_karyotype:
    input:
        seqids="aln/jcvi_karyotype/seqids",
        layouts="aln/jcvi_karyotype/layouts",
        gene_chains="aln/jcvi_karyotype/gene_chains.tsv",
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

        # Composite the ALG colour legend (one brick per shown ALG) onto the
        # right margin, derived from gene_chains.tsv. A no-op if nothing is
        # coloured.
        karyotype_legend \
            --gene-chains gene_chains.tsv \
            --karyotype karyotype.png \
            --layouts layouts
        """


rule jcvi_alg_dotplot:
    """Per-reference dotplot of shared BUSCO genes, colored by ALG identity.

    Positions come from the BUSCO PAF; colors from the karyotype's
    gene_chains.tsv via each line's co:Z: tag; axis order matches the
    karyotype's seqids (off-list chromosomes dropped, for parity with the
    ribbon plot).

    The shared aln_busco.paf has the assembly as query (which dotplotrs draws
    on the Y axis) and each reference as target (X axis). We feed the
    PAF straight through and point the order files at their native axes.

    dotplotrs places the first query (Y) sequence at the top, so we reverse the
    assembly order for the Y axis: the first assembly chromosome sits at the
    bottom, aligning with its matching reference chromosome on the left, giving
    the conventional bottom-left -> top-right diagonal.
    """
    input:
        busco_dir="aln/busco_{accession}",
        gene_chains="aln/jcvi_karyotype/gene_chains.tsv",
        orders="aln/jcvi_karyotype/dotplot_orders",
    output:
        light="aln/jcvi_karyotype/dotplots/{accession}.png",
        dark="aln/jcvi_karyotype/dotplots/{accession}_dark.png",
    benchmark:
        "benchmarks/jcvi_alg_dotplot_{accession}.txt"
    container:
        HOBRAC_TOOLS
    resources:
        mem_mb=8000,
        runtime=30,
    params:
        assembly_key=config["scientific_name"].replace(" ", "_"),
        hide_non_significant=config.get("hide_non_significant", False),
    shell:
        """
        outdir=aln/jcvi_karyotype/dotplots
        mkdir -p $outdir

        # gene -> hex color map (gene_chains.tsv: col3=gene, col4=color).
        id2hex=$(mktemp)
        awk -F'\\t' 'NR>1 && $3!="" {{print $3"\\t"$4}}' {input.gene_chains} > $id2hex

        # Native PAF orientation: assembly = query/Y, reference = target/X.
        ref_order={input.orders}/{wildcards.accession}.order

        # Reverse the assembly order for the Y axis (see rule docstring) so the
        # diagonal runs bottom-left -> top-right.
        assembly_order=$(mktemp)
        awk -F',' '{{for (i = NF; i > 0; i--) printf "%s%s", $i, (i > 1 ? "," : "\\n")}}' \
            {input.orders}/{params.assembly_key}.order > $assembly_order

        hide_flag=""
        if [ "{params.hide_non_significant}" = "True" ]; then
            hide_flag="--hide-unmatched"
        fi

        for theme in light dark; do
            if [ "$theme" = "dark" ]; then
                out={output.dark}
            else
                out={output.light}
            fi
            dotplotrs -p {input.busco_dir}/aln_busco.paf -o $out \
                --colors $id2hex \
                --query-order $assembly_order \
                --target-order $ref_order \
                --theme $theme \
                --line-thickness 10 \
                --font-size 30 \

                $hide_flag
        done

        rm -f $id2hex $assembly_order
        """


rule jcvi_alg_dotplot_grid:
    """Tile every per-reference ALG dotplot into one titled, high-res grid.

    One image per theme, sized to native resolution so the global view stays
    zoomable. Titles resolve from --jcvi-names (full list), else each reference's
    assembly-report organism name, else the accession.
    """
    input:
        dotplots=get_dotplot_grid_inputs,
        order="mash/selected_accessions.txt",
    output:
        light="aln/jcvi_karyotype/dotplots_grid.png",
        dark="aln/jcvi_karyotype/dotplots_grid_dark.png",
    benchmark:
        "benchmarks/jcvi_alg_dotplot_grid.txt"
    container:
        HOBRAC_TOOLS
    resources:
        mem_mb=8000,
        runtime=30,
    params:
        jcvi_names=config.get("jcvi_names", ""),
        assembly_name=config["scientific_name"],
    shell:
        """
        for theme in light dark; do
            if [ "$theme" = "dark" ]; then
                out={output.dark}
            else
                out={output.light}
            fi
            dotplot_grid \
                --dotplots-dir aln/jcvi_karyotype/dotplots \
                --accession-order {input.order} \
                --reference-dir reference \
                --theme $theme \
                -o $out \
                --assembly-name "{params.assembly_name}" \
                --jcvi-names "{params.jcvi_names}"
        done
        """
