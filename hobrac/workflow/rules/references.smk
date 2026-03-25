rule get_reference:
    output: "reference/{accession}.fna"
    container: HOBRAC_TOOLS
    resources:
        mem_mb = 5000,
        runtime = 2 * 60
    benchmark: "benchmarks/get_reference_{accession}.txt"
    shell: """
        cd reference
        find_reference_genomes -d {wildcards.accession} -o {wildcards.accession} 
        mv {wildcards.accession}/*.fna {wildcards.accession}.fna
        rm -r {wildcards.accession}
    """
