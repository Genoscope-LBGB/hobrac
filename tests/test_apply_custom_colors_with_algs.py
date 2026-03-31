from hobrac.jcvi_synteny import (
    BuscoGene,
    PairwiseAssociation,
    apply_custom_colors_with_algs,
)


def _gene(chromosome):
    return BuscoGene(busco_id="", chromosome=chromosome, start=0, end=100)


def test_all_four_cases():
    """Test the four coloring cases"""
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

    chr_to_alg = {
        ("sp1", "chr1"): 0,
        ("sp2", "chrA"): 0,
    }
    alg_colors = {0: "#ff0000"}

    significant_associations = [
        PairwiseAssociation(
            species1="sp1",
            species2="sp2",
            chr1="chr1",
            chr2="chrA",
            p_value=0.001,
            gene_count=10,
        ),
    ]

    custom_colors = {
        "gene_sig_custom": "#00ff00",
        "gene_nonsig_custom": "#0000ff",
    }

    result = apply_custom_colors_with_algs(
        species1_busco,
        species2_busco,
        "sp1",
        "sp2",
        chr_to_alg,
        alg_colors,
        significant_associations,
        custom_colors,
    )

    # Significant pair + in custom file → custom color
    assert result["gene_sig_custom"] == "#00ff00"
    # Significant pair + NOT in custom file → ALG palette color
    assert result["gene_sig_no_custom"] == "#ff0000"
    # Non-significant pair + in custom file → lightgrey
    assert result["gene_nonsig_custom"] == "lightgrey"
    # Non-significant pair + NOT in custom file → lightgrey
    assert result["gene_nonsig_no_custom"] == "lightgrey"


def test_excludes_non_common_genes():
    """Genes not common between both species are excluded."""
    species1_busco = {"common": _gene("chr1"), "only_sp1": _gene("chr1")}
    species2_busco = {"common": _gene("chrA"), "only_sp2": _gene("chrA")}

    result = apply_custom_colors_with_algs(
        species1_busco, species2_busco, "sp1", "sp2", {}, {}, [], {}
    )

    assert "common" in result
    assert "only_sp1" not in result
    assert "only_sp2" not in result


def test_sig_pair_no_alg_mapping_falls_to_lightgrey():
    """Gene on significant pair but chr not in chr_to_alg → lightgrey."""
    species1_busco = {"gene1": _gene("chr1")}
    species2_busco = {"gene1": _gene("chrA")}

    significant_associations = [
        PairwiseAssociation(
            species1="sp1",
            species2="sp2",
            chr1="chr1",
            chr2="chrA",
            p_value=0.001,
            gene_count=10,
        ),
    ]

    result = apply_custom_colors_with_algs(
        species1_busco,
        species2_busco,
        "sp1",
        "sp2",
        {},  # no chr_to_alg mapping
        {},
        significant_associations,
        {},
    )

    assert result["gene1"] == "lightgrey"
