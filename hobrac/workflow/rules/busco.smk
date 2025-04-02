rule get_lineage:
    output: "busco/lineage.txt"
    params:
        taxid = config["taxid"]
    resources:
        threads = 1,
        mem_mb = 5000,
        runtime = 20
    shell: """
        echo {params.taxid} | taxonkit lineage | cut -f 2 > {output}
    """

rule get_busco_datasets:
    input: rules.get_lineage.output
    output: "busco/datasets.txt"
    resources:
        threads = 1,
        mem_mb = 5000,
        runtime = 20
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
        threads = 1,
        mem_mb = 5000,
        runtime = 20
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


rule busco:
    input: 
        rules.get_references.output, 
    output: "busco/busco.done"
    resources:
        threads = 12,
        mem_mb = 60000,
        runtime = 24 * 60
    shell: """
        import glob

        for reference in glob.glob("")
        dataset = 
        busco --metaeuk -i {assembly} -c {cores} -m {mode} --download_path $buscodbpath  -o Busco_{mode}_{busco_db} -l {busco_db}\
    """
                
