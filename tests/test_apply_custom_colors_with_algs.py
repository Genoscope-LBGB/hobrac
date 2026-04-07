from hobrac.jcvi_synteny import (
    BuscoGene,
    apply_custom_colors_with_algs,
)


def _gene(chromosome):
    return BuscoGene(busco_id="", chromosome=chromosome, start=0, end=100)


def test_custom_colors_is_authoritative():
    """With custom_colors, the file drives coloring and chain assignment is ignored.

    Listed genes get their custom color regardless of chain membership; unlisted
    genes always get lightgrey, even if they are on a chain.
    """
    species1_busco = {
        "gene_sig_custom": _gene("chr1"),
        "gene_sig_no_custom": _gene("chr1"),
        "gene_nonsig_custom": _gene("chr3"),
        "gene_nonsig_no_custom": _gene("chr3"),
    }
    species2_busco = {
        "gene_sig_custom": _gene("chrA"),
        "gene_sig_no_custom": _gene("chrA"),
        "gene_nonsig_custom": _gene("chrC"),
        "gene_nonsig_no_custom": _gene("chrC"),
    }

    gene_to_chain = {
        "gene_sig_custom": 0,
        "gene_sig_no_custom": 0,
        "gene_nonsig_custom": -1,
        "gene_nonsig_no_custom": -1,
    }
    chain_colors = {0: "#ff0000"}

    custom_colors = {
        "gene_sig_custom": "#00ff00",
        "gene_nonsig_custom": "#0000ff",
    }

    result = apply_custom_colors_with_algs(
        species1_busco,
        species2_busco,
        gene_to_chain,
        chain_colors,
        custom_colors,
    )

    # In custom file → custom color (chain membership irrelevant)
    assert result["gene_sig_custom"] == "#00ff00"
    assert result["gene_nonsig_custom"] == "#0000ff"
    # NOT in custom file → lightgrey (even when on a chain)
    assert result["gene_sig_no_custom"] == "lightgrey"
    assert result["gene_nonsig_no_custom"] == "lightgrey"


def test_excludes_non_common_genes():
    """Genes not common between both species are excluded."""
    species1_busco = {"common": _gene("chr1"), "only_sp1": _gene("chr1")}
    species2_busco = {"common": _gene("chrA"), "only_sp2": _gene("chrA")}

    result = apply_custom_colors_with_algs(species1_busco, species2_busco, {}, {}, {})

    assert "common" in result
    assert "only_sp1" not in result
    assert "only_sp2" not in result


def test_gene_not_in_gene_to_chain_falls_to_lightgrey():
    """Gene not in gene_to_chain mapping → lightgrey."""
    species1_busco = {"gene1": _gene("chr1")}
    species2_busco = {"gene1": _gene("chrA")}

    result = apply_custom_colors_with_algs(
        species1_busco,
        species2_busco,
        {},  # no gene_to_chain mapping
        {},
        {},
    )

    assert result["gene1"] == "lightgrey"
