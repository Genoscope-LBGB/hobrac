rule get_references:
    output: "references/references.txt"
    params: 
        name = config["scientific_name"]
    resources:
        mem_mb=5000,
        threads=1
    shell: """
        find_reference_genomes -n '{params.name}' -l complete --max-rank order | \
            head -n 6 | \
            tail -n +2 > {output}
    """

rule get_top5_references:
    input: "references/references.txt"
    output: "references/references.ckpt"
    resources:
        mem_mb=5000,
        threads=1
    shell: """
        cd references

        while read line
        do
            name=$(echo $line | cut -f 1 -d ',' | sed 's/ /_/g')
            code=$(echo $line | cut -f 4 -d ',')
            find_reference_genomes -d $code -o $name
        done < references.txt
        
        touch references.ckpt
    """