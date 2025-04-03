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


# checkpoint before_busco:
#     output: touch("busco/collect_references.done")


# def get_references(wildcards):
#     """Function to dynamically list files after checkpoint execution"""
#     files = glob_wildcards("references/{ref_name}.fna").ref_name
#     files = [f.replace(".fna", "") for f in files] 
#     return expand("busco/{ref_name}", zip, ref_name=files) 


# rule launch_buscos:
#     input: get_references


# rule busco:
#     input: 
#         rules.before_busco.output,
#         reference = "references/{ref_name}.fna",
#         dataset = "busco/chosen_dataset.txt"
#     output: directory("busco/{ref_name}")
#     threads: 12
#     resources:
#         mem_mb = 60000,
#         runtime = 24 * 60
#     shell: """
#         dataset=$(cat {input.dataset})

#         cd busco/

#         export buscodbpath="$(pwd)/$(whoami)_buscodb_$$"

#         busco --metaeuk -i ../{input.reference} -c {threads} -m geno \
#             --download_path $buscodbpath  -o {wildcards.ref_name} -l $dataset

#         rm -rf {wildcards.ref_name}/run*/{{busco_sequences,hmmer_output,metaeuk_output,miniprot_output}} ${{buscodbpath}}
#     """
                
