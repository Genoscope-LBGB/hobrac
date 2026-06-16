"""Tests for the dotplot grid: title resolution, organism parsing, grid dims."""

import pytest

from hobrac.jcvi_synteny.grid import (
    grid_dims,
    parse_organism_name,
    resolve_titles,
)

REPORT_HEADER = "\n".join(
    [
        "# Assembly name:  ASM",
        "# Organism name:  Brassica rapa (field mustard)",
        "# Taxid:          3711",
        "# Sequence-Name\tSequence-Role",
    ]
)


def test_parse_organism_name_strips_common_name(tmp_path):
    report = tmp_path / "GCA_1_assembly_report.txt"
    report.write_text(REPORT_HEADER)
    assert parse_organism_name(str(report)) == "Brassica rapa"


def test_parse_organism_name_missing_file_returns_empty(tmp_path):
    assert parse_organism_name(str(tmp_path / "nope.txt")) == ""


def test_parse_organism_name_no_header_returns_empty(tmp_path):
    report = tmp_path / "r.txt"
    report.write_text("# Assembly name:  ASM\n1\tassembled-molecule\n")
    assert parse_organism_name(str(report)) == ""


def test_resolve_titles_full_jcvi_names_used_for_references(tmp_path):
    # 1 assembly + 2 references -> references take names_list[1:].
    titles = resolve_titles(
        ["GCA_1", "GCA_2"],
        jcvi_names="MyAssembly,Species One,Species Two",
        reference_dir=str(tmp_path),
    )
    assert titles == ["Species One", "Species Two"]


def test_resolve_titles_single_name_falls_through_to_report(tmp_path):
    # A single jcvi_name targets the assembly only; references fall through.
    (tmp_path / "GCA_1_assembly_report.txt").write_text(REPORT_HEADER)
    titles = resolve_titles(
        ["GCA_1", "GCA_2"],
        jcvi_names="JustTheAssembly",
        reference_dir=str(tmp_path),
    )
    assert titles == ["Brassica rapa", "GCA_2"]


def test_resolve_titles_no_names_report_then_accession(tmp_path):
    (tmp_path / "GCA_1_assembly_report.txt").write_text(REPORT_HEADER)
    titles = resolve_titles(
        ["GCA_1", "GCA_2"], jcvi_names="", reference_dir=str(tmp_path)
    )
    assert titles == ["Brassica rapa", "GCA_2"]


def test_resolve_titles_mismatched_count_falls_through(tmp_path):
    # Count != 1 + N references -> ignore jcvi_names, use accession fallback.
    titles = resolve_titles(
        ["GCA_1", "GCA_2"], jcvi_names="A,B", reference_dir=str(tmp_path)
    )
    assert titles == ["GCA_1", "GCA_2"]


@pytest.mark.parametrize(
    "n,expected",
    [(1, (1, 1)), (2, (1, 2)), (4, (2, 2)), (5, (2, 3)), (9, (3, 3)), (12, (3, 4))],
)
def test_grid_dims(n, expected):
    assert grid_dims(n) == expected
