rule get_references:
    output: "references/references.txt"
    params: 
        name = config["scientific_name"]
    resources:
        mem_mb = 5000,
        threads = 1,
        runtime = 20
    shell: """
        find_reference_genomes -n '{params.name}' -l complete --max-rank order | \
            head -n 6 | \
            tail -n +2 > {output}
    """

rule get_top5_references:
    input: rules.get_references.output
    output: "references/references.list"
    resources:
        mem_mb = 5000,
        threads = 1,
        runtime = 20
    shell: """
        cd references

        while read line
        do
            name=$(echo $line | cut -f 1 -d ',' | sed 's/ /_/g')
            code=$(echo $line | cut -f 4 -d ',')
            find_reference_genomes -d $code -o $name
        done < references.txt
        
        ls */*.fna > references.list
    """