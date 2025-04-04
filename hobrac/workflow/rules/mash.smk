checkpoint before_mash:
    input: rules.get_top_references.output
    output: touch("mash/collect_references.done")


def get_references():
    import glob, os
    return glob.glob("references/*.fna")


rule launch_mash:
    input: 
        rules.before_mash.output,
        references = get_references(),
        assembly = config["assembly"]
    output: "mash/mash.dist"
    resources:
        mem_mb = 20000,
        runtime = 120
    benchmark: "benchmarks/mash.txt"
    shell: """
        mash sketch -o mash/references {input.references}
        mash info mash/references.msh
        mash dist mash/references.msh {input.assembly} > {output}
    """


rule select_closest_reference:
    input: rules.launch_mash.output
    output: "mash/closest_reference.txt"
    resources:
        mem_mb = 2000,
        runtime = 10
    benchmark: "benchmarks/select_closest_reference.txt"#
    run:
        import os

        with open(input[0]) as mash, open(output[0], "w") as out:
            minimum_distance = 1.0
            path = ""

            for line in mash:
                line = line.rstrip().split("\t")
                distance = float(line[2])
                if distance <= minimum_distance:
                    minimum_distace = distance
                    path = line[0]

            prefix = os.path.basename(path).replace('.fna', '')
            abspath = os.path.abspath(path)
            print(f"{prefix}\t{abspath}", file=out, end="")
            

              
