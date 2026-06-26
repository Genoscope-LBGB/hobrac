import os


rule get_lineage:
    output:
        "busco/lineage.txt",
    benchmark:
        "benchmarks/get_lineage.txt"
    container:
        HOBRAC_TOOLS
    resources:
        mem_mb=5000,
        runtime=20,
    params:
        taxid=config["taxid"],
    shell:
        """
        echo "{params.taxid}" | taxonkit lineage | cut -f 2 > {output}
    """


rule get_busco_datasets:
    input:
        rules.get_lineage.output,
    output:
        "busco/datasets.txt",
    benchmark:
        "benchmarks/get_busco_dataset.txt"
    container:
        "docker://ezlabgva/busco:v6.1.0_cv1"
    resources:
        mem_mb=5000,
        runtime=20,
    shell:
        """
        busco --list-datasets --datasets_version odb12 > {output}
    """


rule get_closest_busco_dataset:
    input:
        lineage=rules.get_lineage.output,
        datasets=rules.get_busco_datasets.output,
    output:
        "busco/chosen_dataset.txt",
    benchmark:
        "benchmarks/get_closest_busco_dataset.txt"
    resources:
        mem_mb=5000,
        runtime=20,
    run:
        import sys

        datasets = {}
        with open(input.datasets[0]) as datasets_file:
            for line in datasets_file:
                line = line.strip()
                line = line.replace(" ", "")
                if not line.startswith("-"):
                    continue
                line = line.replace("-", "")
                dataset_name = line.split("_odb")[0].lower()
                dataset_version = line.lower()
                datasets[dataset_name] = dataset_version.split("[")[0]
        with open(input.lineage[0]) as lineage_file, open(output[0], "w") as out:
            lineage = lineage_file.readline()
            lineage_tree = lineage.split(";")
            # Iterate in reverse order to get the most-specific dataset
            for rank in lineage_tree[::-1]:
                rank = rank.lower()
                if rank in datasets:
                    print(f"{rank}\t{datasets[rank]}", file=out, end="")
                    return
        print(
            "WARNING: no matching dataset found, using eukaryota", file=sys.stderr
        )
        with open(input.lineage[0]) as lineage_file, open(output[0], "w") as out:
            print(f"eukaryota\teukaryota_odb12", file=out, end="")


rule download_busco_dataset:
    input:
        "busco/chosen_dataset.txt",
    output:
        directory("busco/busco_downloads"),
    benchmark:
        "benchmarks/download_busco_dataset.txt"
    container:
        "docker://ezlabgva/busco:v6.1.0_cv1"
    resources:
        mem_mb=5000,
        runtime=60,
    shell:
        """
        version=$(cat {input} | cut -f 2)
        busco --download_path {output} --download $version --datasets_version odb12
    """


rule busco_reference:
    input:
        # ancient(): manual references (-r) are rewritten by main.py on every
        # invocation, bumping the .fna mtime and otherwise re-triggering BUSCO
        # (and a fresh busco_downloads re-download) even when content is
        # unchanged. Ignoring mtime here keeps cached BUSCO results valid. Trade-
        # off: replacing a reference in place no longer auto-reruns BUSCO; delete
        # busco/busco_reference_{accession} to force it.
        fna=ancient("reference/{accession}.fna"),
        dataset="busco/chosen_dataset.txt",
        busco_db="busco/busco_downloads",
    output:
        directory("busco/busco_reference_{accession}"),
    benchmark:
        "benchmarks/busco_reference_{accession}.txt"
    container:
        "docker://ezlabgva/busco:v6.1.0_cv1"
    threads: 12
    resources:
        mem_mb=config["busco_memory"],
        runtime=config["busco_runtime"],
    params:
        method=config["busco_method"],
        fna_path=lambda wildcards, input: (
            input.fna if os.path.isabs(input.fna) else f"../../{input.fna}"
        ),
    shell:
        """
        dataset=$(cat {input.dataset} | cut -f 1)

        # Run in an isolated working directory so concurrent BUSCO jobs don't
        # clobber each other's logs in the shared busco/ folder.
        workdir=busco/tmp_busco_reference_{wildcards.accession}
        rm -rf $workdir
        mkdir -p $workdir
        cd $workdir

        busco --skip_bbtools --{params.method} -i {params.fna_path} -c {threads} -m geno \
            -o busco_reference_{wildcards.accession} -l $dataset \
            --offline --download_path ../busco_downloads --datasets_version odb12

        rm -rf busco_reference_{wildcards.accession}/run*/{{busco_sequences,hmmer_output,metaeuk_output,miniprot_output}}

        cd ..
        rm -rf busco_reference_{wildcards.accession}
        mv tmp_busco_reference_{wildcards.accession}/busco_reference_{wildcards.accession} busco_reference_{wildcards.accession}
        rm -rf tmp_busco_reference_{wildcards.accession}
    """


rule busco_assembly:
    input:
        assembly=config["assembly"],
        dataset=rules.get_closest_busco_dataset.output,
        busco_db="busco/busco_downloads",
    output:
        directory("busco/busco_assembly"),
    benchmark:
        "benchmarks/busco_assembly.txt"
    container:
        "docker://ezlabgva/busco:v6.1.0_cv1"
    threads: 12
    resources:
        mem_mb=config["busco_memory"],
        runtime=config["busco_runtime"],
    params:
        method=config["busco_method"],
        assembly_path=lambda wildcards, input: (
            input.assembly if os.path.isabs(input.assembly) else f"../../{input.assembly}"
        ),
    shell:
        """
        dataset=$(cat {input.dataset} | cut -f 1)

        # Run in an isolated working directory so concurrent BUSCO jobs don't
        # clobber each other's logs in the shared busco/ folder.
        workdir=busco/tmp_busco_assembly
        rm -rf $workdir
        mkdir -p $workdir
        cd $workdir

        busco --skip_bbtools --{params.method} -i {params.assembly_path} -c {threads} -m geno \
            -o busco_assembly -l $dataset \
            --offline --download_path ../busco_downloads --datasets_version odb12

        rm -rf busco_assembly/run*/{{busco_sequences,hmmer_output,metaeuk_output,miniprot_output}}

        cd ..
        rm -rf busco_assembly
        mv tmp_busco_assembly/busco_assembly busco_assembly
        rm -rf tmp_busco_assembly
    """


rule cleanup_busco_downloads:
    input:
        get_pipeline_targets,
    output:
        touch("aln/cleanup.done"),
    resources:
        mem_mb=1000,
        runtime=5,
    shell:
        """
        rm -rf busco/busco_downloads reference/*.fna assembly/*.fna
    """


rule busco_to_paf:
    input:
        reference=ancient("reference/{accession}.fna"),
        assembly=config["assembly"],
        busco_reference=lambda wildcards: config.get(
            "busco_reference_override", f"busco/busco_reference_{wildcards.accession}"
        ),
        busco_assembly=lambda wildcards: config.get(
            "busco_assembly_override", "busco/busco_assembly"
        ),
    output:
        directory("aln/busco_{accession}"),
    benchmark:
        "benchmarks/busco_to_paf_{accession}.txt"
    container:
        HOBRAC_TOOLS
    resources:
        mem_mb=50000,
        runtime=600,
    params:
        prefix_assembly=config["scientific_name"].replace(" ", "_"),
    shell:
        """
        busco_to_paf --busco_query {input.busco_assembly}/run*/full_table.tsv \
            --busco_ref {input.busco_reference}/run*/full_table.tsv \
            --query {input.assembly} --ref {input.reference} --out {output}

        mv {output}/query_assembly.idx "{output}/busco_query_{params.prefix_assembly}.idx"
        mv {output}/target_reference.idx {output}/busco_target_{wildcards.accession}.idx

        dotplotrs -p {output}/aln_busco.paf -o {output}/dotplot_busco.png --line-thickness 4 --dump-significance {output}/significance_results.txt
        dotplotrs -p {output}/aln_busco.paf -o {output}/dotplot_busco_bw.png --line-thickness 4 --no-color
    """
