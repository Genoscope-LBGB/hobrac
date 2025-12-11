rule download_db:
    output: "mash/mash_db.msh"
    params: taxid = config["taxid"]
    container: "docker://ghcr.io/cea-lbgb/hobrac-tools:latest"
    benchmark: "benchmarks/download_db.txt"
    resources:
        mem_mb = 10000,
        runtime = 30
    shell: """
        phylum=$(echo {params.taxid} | taxonkit reformat -I 1 --format '{{p}}' -r 'no_returned_phylum' | cut -f 2)
        wget https://www.genoscope.cns.fr/lbgb/mash/${{phylum}}.msh -O {output}
        echo ${{phylum}} > {output}.phylum
    """


rule launch_mash:
    input: 
        mashdb = rules.download_db.output,
        assembly = config["assembly"]
    output: "mash/mash.dist"
    container: "docker://ghcr.io/cea-lbgb/hobrac-tools:latest"
    resources:
        mem_mb = 20000,
        runtime = 120
    benchmark: "benchmarks/mash.txt"
    shell: """
        mash info {input.mashdb} > {input.mashdb}.info
        mash dist -s 10000 {input.mashdb} {input.assembly} > {output}
    """


rule select_closest_reference:
    input: rules.launch_mash.output
    output: "mash/closest_reference.txt"
    params:
        allow_zero_distance = config["allow_zero_distance"],
        allow_same_taxid = config["allow_same_taxid"],
        taxid = config["taxid"]
    resources:
        mem_mb = 2000,
        runtime = 10
    benchmark: "benchmarks/select_closest_reference.txt"
    run:
        import os

        with open(input[0]) as mash, open(output[0], "w") as out:
            minimum_distance = 1.0
            accession = ""

            for line in mash:
                line = line.rstrip().split("\t")
                distance = float(line[2])

                taxid = "0"
                if ":" in line[0]:
                    taxid = line[0].split(":")[1]
                if not params.allow_same_taxid and params.taxid == taxid:
                    continue

                if distance == 0.0 and not params.allow_zero_distance:
                    continue
                    
                if distance <= minimum_distance:
                    minimum_distance = distance
                    accession = line[0].split(":")[0]

            print(accession, file=out, end="")
