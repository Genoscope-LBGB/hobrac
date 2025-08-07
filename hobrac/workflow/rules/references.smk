rule get_top_reference:
    input: rules.select_closest_reference.output
    output: "reference/reference.txt"
    log: "logs/references/get_top_reference.log"
    resources:
        mem_mb = 5000,
        runtime = 2 * 60
    benchmark: "benchmarks/get_top_reference.txt"
    shell: """
        accession=$(cat {input})
        
        cd reference 2>> {log}
        find_reference_genomes -d $accession -o $accession 2>> {log}
        mv $accession/*.fna $accession.fna 2>> {log}
        rm -r $accession 2>> {log}
        
        ls $(pwd)/*.fna > reference.txt
    """
