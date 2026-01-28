import os


rule get_lineage:
    output: "busco/lineage.txt"
    params:
        taxid = config["taxid"]
    container: HOBRAC_TOOLS
    resources:
        mem_mb = 5000,
        runtime = 20
    benchmark: "benchmarks/get_lineage.txt"
    shell: """
        echo {params.taxid} | taxonkit lineage | cut -f 2 > {output}
    """


rule get_busco_datasets:
    input: rules.get_lineage.output
    output: "busco/datasets.txt"
    container: "docker://ezlabgva/busco:v6.0.0_cv1"
    resources:
        mem_mb = 5000,
        runtime = 20
    benchmark: "benchmarks/get_busco_dataset.txt"
    shell: """
        busco --list-datasets > {output}
    """


rule get_closest_busco_dataset:
    input:
        lineage = rules.get_lineage.output,
        datasets = rules.get_busco_datasets.output
    output: "busco/chosen_dataset.txt"
    resources:
        mem_mb = 5000,
        runtime = 20
    benchmark: "benchmarks/get_closest_busco_dataset.txt"
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

        print("WARNING: no matching dataset found, using eukaryota", file=sys.stderr)
        with open(input.lineage[0]) as lineage_file, open(output[0], "w") as out:
            print(f"eukaryota\teukaryota_odb12", file=out, end="")


rule download_busco_dataset:
    input: "busco/chosen_dataset.txt"
    output: directory("busco/busco_downloads")
    container: "docker://ezlabgva/busco:v6.0.0_cv1"
    resources:
        mem_mb = 5000,
        runtime = 60
    benchmark: "benchmarks/download_busco_dataset.txt"
    shell: """
        version=$(cat {input} | cut -f 2)
        busco --download_path {output} --download $version
    """


rule busco_reference:
    input:
        fna = "reference/{accession}.fna",
        dataset = "busco/chosen_dataset.txt"
    output: directory("busco/busco_reference_{accession}")
    params:
        method = config["busco_method"],
        fna_path = lambda wildcards, input: input.fna if os.path.isabs(input.fna) else f"../{input.fna}"
    threads: 12
    container: "docker://ezlabgva/busco:v6.0.0_cv1"
    resources:
        mem_mb = config["busco_memory"],
        runtime = config["busco_runtime"]
    benchmark: "benchmarks/busco_reference_{accession}.txt"
    shell: """
        dataset=$(cat {input.dataset} | cut -f 1)

        cd busco/

        busco --skip_bbtools --{params.method} -i {params.fna_path} -c {threads} -m geno \
            -o busco_reference_{wildcards.accession} -l $dataset \
            --offline --download_path ../busco_downloads

        rm -rf busco_reference_{wildcards.accession}/run*/{{busco_sequences,hmmer_output,metaeuk_output,miniprot_output}}
    """
    
    
rule busco_assembly:
    input:
        assembly = config["assembly"],
        dataset = rules.get_closest_busco_dataset.output
    output: directory("busco/busco_assembly")
    params:
        method = config["busco_method"],
        assembly_path = lambda wildcards, input: input.assembly if os.path.isabs(input.assembly) else f"../{input.assembly}"
    threads: 12
    container: "docker://ezlabgva/busco:v6.0.0_cv1"
    resources:
        mem_mb = config["busco_memory"],
        runtime = config["busco_runtime"]
    benchmark: "benchmarks/busco_assembly.txt"
    shell: """
        dataset=$(cat {input.dataset} | cut -f 1)

        cd busco/

        busco --skip_bbtools --{params.method} -i {params.assembly_path} -c {threads} -m geno \
            -o busco_assembly -l $dataset \
            --offline --download_path ../busco_downloads

        rm -rf busco_assembly/run*/{{busco_sequences,hmmer_output,metaeuk_output,miniprot_output}}
    """


rule busco_to_paf:
    input:
        reference = "reference/{accession}.fna",
        assembly = config["assembly"],
        busco_reference = lambda wildcards: config.get("busco_reference_override", f"busco/busco_reference_{wildcards.accession}"),
        busco_assembly = lambda wildcards: config.get("busco_assembly_override", "busco/busco_assembly")
    output: directory("aln/busco_{accession}")
    params:
        prefix_assembly = config["scientific_name"].replace(" ", "_")
    container: HOBRAC_TOOLS
    resources:
        mem_mb = 50000,
        runtime = 600
    benchmark: "benchmarks/busco_to_paf_{accession}.txt"
    shell: """
        busco_to_paf --busco_query {input.busco_assembly}/run*/full_table.tsv \
            --busco_ref {input.busco_reference}/run*/full_table.tsv \
            --query {input.assembly} --ref {input.reference} --out {output}

        mv {output}/query_assembly.idx {output}/busco_query_{params.prefix_assembly}.idx
        mv {output}/target_reference.idx {output}/busco_target_{wildcards.accession}.idx

        dotplotrs -p {output}/aln_busco.paf -o {output}/busco_significance.png --line-thickness 8
        dotplotrs -p {output}/aln_busco.paf -o {output}/busco_gravity.png --line-thickness 8 --gravity-ordering-only
        dotplotrs -p {output}/aln_busco.paf -o {output}/busco_bw.png --line-thickness 8 --no-color
        dotplotrs -p {output}/aln_busco.paf -o {output}/busco_bw_gravity.png --line-thickness 8 --no-color --gravity-ordering-only
    """
