import os


rule get_lineage:
    output: "busco/lineage.txt"
    params:
        taxid = config["taxid"]
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
    resources:
        mem_mb = 5000,
        runtime = 20
    benchmark: "benchmarks/get_busco_dataset.txt"
    shell: """
        busco --list-datasets > {output}
        rm -rf busco_downloads
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

        datasets = set()
        with open(input.datasets[0]) as datasets_file:
            for line in datasets_file:
                line = line.strip()
                line = line.replace(" ", "")

                if not line.startswith("-"):
                    continue

                line = line.replace("-", "")
                dataset_name = line.split("_odb")[0].lower()
                datasets.add(dataset_name)

        with open(input.lineage[0]) as lineage_file, open(output[0], "w") as out:
            lineage = lineage_file.readline()
            lineage_tree = lineage.split(";")

            # Iterate in reverse order to get the most-specific dataset
            for rank in lineage_tree[::-1]:
                rank = rank.lower()
                if rank in datasets:
                    print(rank, file=out, end="")
                    return

        print("ERROR: no matching dataset found", file=sys.stderr)
        sys.exit(1)     


rule busco_reference:
    input: 
        mash_output = rules.select_closest_reference.output,
        dataset = "busco/chosen_dataset.txt"
    output: directory("busco/busco_reference")
    threads: 12
    resources:
        mem_mb = config["busco_memory"],
        runtime = 24 * 60
    shell: """
        dataset=$(cat {input.dataset})
        prefix=$(cat {input.mash_output} | cut -f 1)
        filepath=$(cat {input.mash_output} | cut -f 2)

        cd busco/

        export buscodbpath="$(pwd)/$(whoami)_buscodb_$$"

        busco --metaeuk -i $filepath -c {threads} -m geno \
            --download_path $buscodbpath  -o busco_reference -l $dataset
        
        ln -s busco_reference busco_$prefix
        rm -rf busco_reference/run*/{{busco_sequences,hmmer_output,metaeuk_output,miniprot_output}} ${{buscodbpath}}
    """


rule busco_assembly:
    input:
        assembly = config["assembly"],
        dataset = rules.get_closest_busco_dataset.output
    output: directory("busco/busco_assembly")
    threads: 12
    resources:
        mem_mb = config["busco_memory"],
        runtime = 24 * 60
    shell: """
        dataset=$(cat {input.dataset})

        cd busco/

        export buscodbpath="$(pwd)/$(whoami)_buscodb_$$"

        busco --metaeuk -i {input.assembly} -c {threads} -m geno \
            --download_path $buscodbpath  -o busco_assembly -l $dataset
        
        rm -rf busco_assembly/run*/{{busco_sequences,hmmer_output,metaeuk_output,miniprot_output}} ${{buscodbpath}}
    """
                

rule busco_to_paf:
    input:
        mash_output = rules.select_closest_reference.output,
        assembly = config["assembly"],
        busco_reference = "busco/busco_reference",
        busco_assembly = "busco/busco_assembly"
    output: directory("aln/busco")
    params:
        prefix_assembly = config["scientific_name"].replace(" ", "_")
    resources:
        mem_mb = 10000,
        runtime = 20
    benchmark: "benchmarks/busco_to_paf.txt"
    shell: """
        prefix_ref=$(cat {input.mash_output} | cut -f 1)
        filepath=$(cat {input.mash_output} | cut -f 2)

        busco_to_paf --busco_query {input.busco_assembly}/run*/full_table.tsv \
            --busco_ref {input.busco_reference}/run*/full_table.tsv \
            --query {input.assembly} --ref $filepath --out {output}

        mv {output}/query_assembly.idx {output}/busco_query_{params.prefix_assembly}.idx
        mv {output}/target_reference.idx {output}/busco_target_${{prefix_ref}}.idx

        dotplotrs -p {output}/aln_busco.paf -o {output}/busco.png --line-thickness 4
    """