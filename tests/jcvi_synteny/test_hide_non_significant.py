"""Tests verifying --jcvi-hide-non-significant works across all flag combinations."""

import os
import tempfile

import pytest

from hobrac.jcvi_synteny.coloring import (
    apply_custom_colors,
    apply_custom_colors_with_algs,
)
from hobrac.jcvi_synteny.models import BuscoGene
from hobrac.jcvi_synteny.output import generate_links_file


def _gene(chromosome, start):
    """Create a BuscoGene with distinct positions for block detection."""
    return BuscoGene(busco_id="", chromosome=chromosome, start=start, end=start + 100)


@pytest.fixture
def species_busco():
    """Each gene on a unique chromosome pair to avoid block merging."""
    sp1 = {
        "sig_gene1": _gene("chr1", 1000),
        "sig_gene2": _gene("chr2", 1000),
        "nonsig_gene1": _gene("chr3", 1000),
        "nonsig_gene2": _gene("chr4", 1000),
    }
    sp2 = {
        "sig_gene1": _gene("chrA", 1000),
        "sig_gene2": _gene("chrB", 1000),
        "nonsig_gene1": _gene("chrC", 1000),
        "nonsig_gene2": _gene("chrD", 1000),
    }
    return sp1, sp2


@pytest.fixture
def chain_data(species_busco):
    sp1, sp2 = species_busco
    gene_to_chain = {
        "sig_gene1": 0,
        "sig_gene2": 0,
        "nonsig_gene1": -1,
        "nonsig_gene2": -1,
    }
    chain_colors = {0: "#ff0000"}
    return gene_to_chain, chain_colors


@pytest.fixture
def custom_colors():
    return {"sig_gene1": "#00ff00", "nonsig_gene1": "#0000ff"}


def _read_simple_file(path):
    """Read .simple file and return list of (start_gene, color) tuples."""
    lines = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            gene1_start = parts[0]
            # Parse color from "color*sp_gene" format
            if "*" in gene1_start:
                color, gene = gene1_start.split("*", 1)
            else:
                color, gene = "lightgrey", gene1_start
            lines.append((gene, color))
    return lines


def _write_links(sp1_busco, sp2_busco, gene_colors, hide):
    """Write links to a temp file and return parsed lines."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".simple", delete=False) as f:
        path = f.name
    try:
        generate_links_file(
            sp1_busco,
            sp2_busco,
            gene_colors,
            "sp1",
            "sp2",
            path,
            hide_non_significant=hide,
        )
        return _read_simple_file(path)
    finally:
        os.unlink(path)


# Row 1: No custom_colors, hide=False → all links shown, non-sig in lightgrey
def test_no_custom_hide_false(species_busco, chain_data):
    sp1, sp2 = species_busco
    gene_to_chain, chain_colors = chain_data

    gene_colors = apply_custom_colors_with_algs(
        sp1, sp2, gene_to_chain, chain_colors, {}
    )
    lines = _write_links(sp1, sp2, gene_colors, hide=False)

    sig_blocks = [g for g, c in lines if c != "lightgrey"]
    nonsig_blocks = [g for g, c in lines if c == "lightgrey"]
    assert len(sig_blocks) > 0
    assert len(nonsig_blocks) > 0


# Row 2: No custom_colors, hide=True → only ALG-significant links shown
def test_no_custom_hide_true(species_busco, chain_data):
    sp1, sp2 = species_busco
    gene_to_chain, chain_colors = chain_data

    gene_colors = apply_custom_colors_with_algs(
        sp1, sp2, gene_to_chain, chain_colors, {}
    )
    lines = _write_links(sp1, sp2, gene_colors, hide=True)

    # Only significant (non-lightgrey) blocks remain
    assert len(lines) > 0
    assert all(color != "lightgrey" for _, color in lines)


# Row 3: custom_colors + skip_alg=False + hide=False → all links, sig use custom/ALG
def test_custom_with_alg_hide_false(species_busco, chain_data, custom_colors):
    sp1, sp2 = species_busco
    gene_to_chain, chain_colors = chain_data

    gene_colors = apply_custom_colors_with_algs(
        sp1, sp2, gene_to_chain, chain_colors, custom_colors
    )
    lines = _write_links(sp1, sp2, gene_colors, hide=False)

    # Both colored and grey blocks present
    non_grey = [c for _, c in lines if c != "lightgrey"]
    grey = [c for _, c in lines if c == "lightgrey"]
    assert len(non_grey) > 0
    assert len(grey) > 0


# Row 4: custom_colors + skip_alg=False + hide=True → only significant links
def test_custom_with_alg_hide_true(species_busco, chain_data, custom_colors):
    sp1, sp2 = species_busco
    gene_to_chain, chain_colors = chain_data

    gene_colors = apply_custom_colors_with_algs(
        sp1, sp2, gene_to_chain, chain_colors, custom_colors
    )
    lines = _write_links(sp1, sp2, gene_colors, hide=True)

    # Only non-lightgrey blocks remain
    assert len(lines) > 0
    assert all(color != "lightgrey" for _, color in lines)
    # Custom color for sig_gene1 should be present
    assert any("#00ff00" in color for _, color in lines)


# Row 5: custom_colors + skip_alg=True + hide=False → all links, unlisted grey
def test_custom_skip_alg_hide_false(species_busco, custom_colors):
    sp1, sp2 = species_busco

    gene_colors = apply_custom_colors(sp1, sp2, custom_colors)
    lines = _write_links(sp1, sp2, gene_colors, hide=False)

    # Both custom-colored and grey blocks present
    non_grey = [c for _, c in lines if c != "lightgrey"]
    grey = [c for _, c in lines if c == "lightgrey"]
    assert len(non_grey) > 0
    assert len(grey) > 0


# Row 6: custom_colors + skip_alg=True + hide=True → only custom-colored links
def test_custom_skip_alg_hide_true(species_busco, custom_colors):
    sp1, sp2 = species_busco

    gene_colors = apply_custom_colors(sp1, sp2, custom_colors)
    lines = _write_links(sp1, sp2, gene_colors, hide=True)

    # Only custom-colored (non-lightgrey) blocks remain
    assert len(lines) > 0
    assert all(color != "lightgrey" for _, color in lines)


# custom + ALG + hide → only genes that are BOTH on a chain AND listed in the
# custom file survive. Everything else (off-chain or unlisted) is lightgrey
# and therefore filtered out.
def test_only_on_chain_and_listed_survives_when_hidden(
    species_busco, chain_data, custom_colors
):
    sp1, sp2 = species_busco
    gene_to_chain, chain_colors = chain_data

    gene_colors = apply_custom_colors_with_algs(
        sp1, sp2, gene_to_chain, chain_colors, custom_colors
    )
    lines = _write_links(sp1, sp2, gene_colors, hide=True)

    all_colors = {color for _, color in lines}
    # sig_gene1 is on chain AND listed → survives with custom color.
    assert "#00ff00" in all_colors
    # nonsig_gene1 is listed but off-chain → lightgrey → hidden.
    assert "#0000ff" not in all_colors
    # sig_gene2 is on chain but unlisted → lightgrey → hidden.
    # (The chain palette color must not appear either.)
    assert "#ff0000" not in all_colors
