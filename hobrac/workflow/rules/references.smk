rule get_references:
    output: "references/references.txt"
    params: 
        name = config["scientific_name"],
        taxid = config["taxid"],
        allow_same_taxid = config["allow_same_taxid"],
        nblines = 5 + 1
    resources:
        mem_mb = 5000,
        runtime = 20
    benchmark: "benchmarks/get_references.txt"
    shell: """
        find_reference_genomes -n '{params.name}' -l complete --max-rank phylum > genomes.txt

        if [[ '{params.allow_same_taxid}' == 'False' ]]; then
            cat genomes.txt | \
                awk -v FS=',' 'BEGIN{{nb=0}} {{ if($2 != {params.taxid} && nb < {params.nblines}) {{ print $0; nb+=1}} }}' | \
                tail -n +2 > {output}
        else
            cat genomes.txt | \
                    awk -v FS=',' 'BEGIN{{nb=0}} {{ if(nb < {params.nblines}) {{ print $0; nb+=1}} }}' | \
                    tail -n +2 > {output}
        fi

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
            code=$(echo $line | cut -f 4 -d ',')
            find_reference_genomes -d $code -o $name 2> find_reference_genomes.err
            mv $name/*.fna $name.fna
            rm -r $name
        done < references.txt
        
        ls *.fna > references.list
    """
