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
        
        cd reference
        find_reference_genomes -d $accession -o $accession 
        mv $accession/*.fna $accession.fna
        rm -r $accession
        
        ls $(pwd)/*.fna > reference.txt
    """
