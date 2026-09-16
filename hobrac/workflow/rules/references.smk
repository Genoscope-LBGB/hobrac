rule get_reference:
    output:
        fna="reference/{accession}.fna",
        report="reference/{accession}_assembly_report.txt",
    benchmark:
        "benchmarks/get_reference_{accession}.txt"
    container:
        HOBRAC_TOOLS
    resources:
        mem_mb=5000,
        runtime=2 * 60,
    params:
        rename_flag="" if config.get("busco_reference_override") else "-r",
    shell:
        """
        cd reference
        find_reference_genomes -d {wildcards.accession} -o {wildcards.accession} {params.rename_flag}
        mv {wildcards.accession}/*.fna {wildcards.accession}.fna
        mv {wildcards.accession}/*_assembly_report.txt {wildcards.accession}_assembly_report.txt
        rm -r {wildcards.accession}
    """
