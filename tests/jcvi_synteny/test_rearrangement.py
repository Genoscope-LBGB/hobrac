import pytest

from hobrac.jcvi_synteny.io import parse_custom_algs
from hobrac.jcvi_synteny.models import BuscoGene
from hobrac.jcvi_synteny.rearrangement import (
    calculate_rearrangement_indices,
    summarize_rearrangement_indices,
)


def _gene(chromosome, start=0):
    return BuscoGene(busco_id="", chromosome=chromosome, start=start, end=start + 100)


def test_no_rearrangement_index_is_zero():
    busco = {
        "a1": _gene("chr1"),
        "a2": _gene("chr1"),
        "b1": _gene("chr2"),
        "b2": _gene("chr2"),
    }
    gene_to_alg = {"a1": "ALG_A", "a2": "ALG_A", "b1": "ALG_B", "b2": "ALG_B"}

    rows = calculate_rearrangement_indices("sp", busco, gene_to_alg)
    ri, splitting_index, combining_index, n_algs, n_genes = (
        summarize_rearrangement_indices(rows)
    )

    assert ri == 0
    assert splitting_index == 0
    assert combining_index == 0
    assert n_algs == 2
    assert n_genes == 4


def test_combining_by_chromosome_fusion_matches_paper_example():
    busco = {
        "a1": _gene("chr1"),
        "a2": _gene("chr1"),
        "b1": _gene("chr1"),
        "b2": _gene("chr1"),
    }
    gene_to_alg = {"a1": "ALG_A", "a2": "ALG_A", "b1": "ALG_B", "b2": "ALG_B"}

    rows = calculate_rearrangement_indices("sp", busco, gene_to_alg)
    ri, splitting_index, combining_index, _, _ = summarize_rearrangement_indices(rows)

    assert ri == 0.5
    assert splitting_index == 0
    assert combining_index == 0.5


def test_splitting_and_combining_use_best_chromosome():
    busco = {
        "a1": _gene("chr1"),
        "a2": _gene("chr1"),
        "a3": _gene("chr2"),
        "b1": _gene("chr1"),
    }
    gene_to_alg = {
        "a1": "ALG_A",
        "a2": "ALG_A",
        "a3": "ALG_A",
        "b1": "ALG_B",
    }

    rows = calculate_rearrangement_indices("sp", busco, gene_to_alg)
    by_alg = {row.alg: row for row in rows}

    alg_a = by_alg["ALG_A"]
    assert alg_a.best_chromosome == "chr1"
    assert alg_a.splitting_parameter == 2 / 3
    assert alg_a.combining_parameter == 2 / 3
    assert alg_a.rearrangement_index == 1 - ((2 / 3) * (2 / 3))


def test_unannotated_and_absent_algs_are_not_scored():
    busco = {
        "a1": _gene("chr1"),
        "unannotated": _gene("chr2"),
    }
    gene_to_alg = {
        "a1": "ALG_A",
        "missing_from_species": "ALG_B",
    }

    rows = calculate_rearrangement_indices("sp", busco, gene_to_alg)

    assert [row.alg for row in rows] == ["ALG_A"]
    assert summarize_rearrangement_indices(rows)[0] == 0


def test_parse_custom_algs_uses_third_column(tmp_path):
    color_file = tmp_path / "colors.tsv"
    color_file.write_text(
        "# BUSCO_ID color ALG\n"
        "g1\t141,78,106\tALG_A\n"
        "g2\t1,2,3\tALG_B\n"
    )

    assert parse_custom_algs(str(color_file)) == {
        "g1": "ALG_A",
        "g2": "ALG_B",
    }


def test_parse_custom_algs_warns_on_wrong_column_count(tmp_path):
    color_file = tmp_path / "colors.tsv"
    color_file.write_text(
        "g1\t141,78,106\tALG_A\n"
        "g2\t1,2,3\n"
        "g3\t4,5,6\tALG_C\textra\n"
    )

    with pytest.warns(UserWarning, match="expected 3 columns") as warnings:
        parsed = parse_custom_algs(str(color_file))

    assert len(warnings) == 2
    assert parsed == {
        "g1": "ALG_A",
        "g3": "ALG_C",
    }
