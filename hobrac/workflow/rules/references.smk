rule get_references:
    output: "references/references.txt"
    params: 
        name = config["scientific_name"],
        nblines = 5 + 1
    resources:
        mem_mb = 5000,
        runtime = 20
    benchmark: "benchmarks/get_references.txt"
    shell: """
        find_reference_genomes -n '{params.name}' -l complete --max-rank class > genomes.txt

        cat genomes.txt | \
            head -n {params.nblines} | \
            tail -n +2 > {output}

        rm genomes.txt
    """


rule get_top_references:
    input: rules.get_references.output
    output: "references/references.list"
    resources:
        mem_mb = 5000,
        runtime = 2 * 60
    benchmark: "benchmarks/get_top_references.txt"
    shell: """
        cd references

        while read line
        do
            name=$(echo $line | cut -f 1 -d ',' | sed 's/ /_/g')
            code=$(echo $line | cut -f 3 -d ',')
            find_reference_genomes -d $code -o $name 2> find_reference_genomes.err
            mv $name/*.fna $name.fna
            rm -r $name
        done < references.txt
        
        ls *.fna > references.list
    """