import csv
import math

import pytest

from hobrac.jcvi_synteny.io import parse_custom_algs
from hobrac.jcvi_synteny.models import BuscoGene
from hobrac.jcvi_synteny.rearrangement import (
    calculate_rearrangement_indices,
    save_rearrangement_indices,
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
        summarize_rearrangement_indices(rows, min_ri_alg_genes=1)
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
    ri, splitting_index, combining_index, _, _ = summarize_rearrangement_indices(
        rows, min_ri_alg_genes=1
    )

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
    assert summarize_rearrangement_indices(rows, min_ri_alg_genes=1)[0] == 0


def test_parse_custom_algs_uses_third_column(tmp_path):
    color_file = tmp_path / "colors.tsv"
    color_file.write_text(
        "# BUSCO_ID color ALG\ng1\t141,78,106\tALG_A\ng2\t1,2,3\tALG_B\n"
    )

    assert parse_custom_algs(str(color_file)) == {
        "g1": "ALG_A",
        "g2": "ALG_B",
    }


def test_parse_custom_algs_warns_on_wrong_column_count(tmp_path):
    color_file = tmp_path / "colors.tsv"
    color_file.write_text("g1\t141,78,106\tALG_A\ng2\t1,2,3\ng3\t4,5,6\tALG_C\textra\n")

    with pytest.warns(UserWarning, match="expected 3 columns") as warnings:
        parsed = parse_custom_algs(str(color_file))

    assert len(warnings) == 2
    assert parsed == {
        "g1": "ALG_A",
        "g3": "ALG_C",
    }


@pytest.fixture
def supported_rows():
    busco = {}
    annotations = {}
    for alg, count, chrom in [("A", 5, "chr1"), ("B", 4, "chr1"), ("C", 6, "chr2")]:
        for i in range(count):
            gene_id = f"{alg}{i}"
            busco[gene_id] = _gene(chrom)
            annotations[gene_id] = alg
    return calculate_rearrangement_indices("sp", busco, annotations)


def test_default_threshold_is_inclusive_and_preserves_denominators(supported_rows):
    ri, splitting, combining, n_algs, n_genes = summarize_rearrangement_indices(
        iter(supported_rows)
    )
    # A qualifies at exactly five genes; B's four genes still count on chr1.
    assert ri == pytest.approx((1 - 5 / 9) / 2)
    assert combining == pytest.approx(ri)
    assert splitting == 0
    assert (n_algs, n_genes) == (2, 11)
    assert supported_rows[0].chromosome_gene_count == 9


def test_custom_threshold_and_legacy_threshold(supported_rows):
    assert summarize_rearrangement_indices(supported_rows, 6) == (0, 0, 0, 1, 6)
    ri, _, _, n_algs, n_genes = summarize_rearrangement_indices(supported_rows, 1)
    assert ri == pytest.approx(
        sum(row.rearrangement_index for row in supported_rows) / 3
    )
    assert (n_algs, n_genes) == (3, 15)


@pytest.mark.parametrize("empty_input", [False, True])
def test_no_eligible_algs_have_undefined_indices(supported_rows, empty_input):
    rows = [] if empty_input else supported_rows
    ri, splitting, combining, n_algs, n_genes = summarize_rearrangement_indices(rows, 7)
    assert all(math.isnan(value) for value in (ri, splitting, combining))
    assert (n_algs, n_genes) == (0, 0)


@pytest.mark.parametrize("threshold", [0, -1])
def test_invalid_threshold_rejected(threshold):
    with pytest.raises(ValueError, match="must be at least 1"):
        summarize_rearrangement_indices([], threshold)


def _read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_saved_scores_and_warnings_are_species_specific(tmp_path, supported_rows):
    summary = tmp_path / "summary.tsv"
    detail = tmp_path / "detail.tsv"
    warnings = tmp_path / "warnings.tsv"
    save_rearrangement_indices(
        {"sp": supported_rows, "empty": []},
        str(summary),
        str(detail),
        str(warnings),
        expected_algs=["A", "B", "C", "D", "D"],
    )
    species = {row["species"]: row for row in _read_tsv(summary)}
    assert species["sp"]["n_algs"] == "2"
    assert species["sp"]["n_genes"] == "11"
    assert float(species["sp"]["rearrangement_index"]) == pytest.approx(2 / 9, abs=1e-6)
    assert species["empty"]["rearrangement_index"] == "nan"
    assert species["empty"]["n_algs"] == "0"
    details = _read_tsv(detail)
    assert [(row["alg"], row["included_in_ri"]) for row in details] == [
        ("A", "true"),
        ("B", "false"),
        ("C", "true"),
    ]
    assert details[0]["chromosome_gene_count"] == "9"
    warning_rows = _read_tsv(warnings)
    assert [
        (row["species"], row["alg"], row["alg_gene_count"]) for row in warning_rows
    ] == [
        ("empty", "A", "0"),
        ("empty", "B", "0"),
        ("empty", "C", "0"),
        ("empty", "D", "0"),
        ("sp", "B", "4"),
        ("sp", "D", "0"),
    ]
    assert all(row["min_ri_alg_genes"] == "5" for row in warning_rows)
    assert all("Excluded from species RI" in row["warning"] for row in warning_rows)

    # Reruns replace obsolete warnings; no exclusions leave a header-only file.
    save_rearrangement_indices(
        {"sp": supported_rows},
        str(summary),
        str(detail),
        str(warnings),
        min_ri_alg_genes=1,
        expected_algs=["A", "B", "C"],
    )
    assert _read_tsv(warnings) == []
    assert warnings.read_text().startswith("species\talg\talg_gene_count\t")
    assert all(row["included_in_ri"] == "true" for row in _read_tsv(detail))
