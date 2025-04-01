rule get_lineage:
    output: "busco/lineage.txt"
    params:
        taxid = config["taxid"]
    resources:
        threads = 1,
        mem_mb = 5000
    shell:"""
        echo {params.taxid} | taxonkit lineage | cut -f 2 > {output}
    """