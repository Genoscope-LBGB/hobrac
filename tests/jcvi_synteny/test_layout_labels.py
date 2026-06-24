"""Tests for karyotype track labels (output.resolve_layout_labels).

These mirror tests/jcvi_synteny/test_dotplot_grid.py: the karyotype labels must
resolve the same way the dotplot grid titles do so the two figures agree.
"""

from hobrac.jcvi_synteny.output import resolve_layout_labels

REPORT_HEADER = "\n".join(
    [
        "# Assembly name:  ASM",
        "# Organism name:  Brassica rapa (field mustard)",
        "# Taxid:          3711",
        "# Sequence-Name\tSequence-Role",
    ]
)


def _keys():
    """Stable keys: underscored assembly name, then reference accessions."""
    return ["Felimare_picta", "GCA_1", "GCA_2"]


def test_no_names_uses_display_name_then_organism_then_accession(tmp_path):
    (tmp_path / "GCA_1_assembly_report.txt").write_text(REPORT_HEADER)
    labels = resolve_layout_labels(
        _keys(),
        assembly_display_name="Felimare picta",
        custom_names="",
        reference_dir=str(tmp_path),
    )
    # Assembly: spaced display name (not the underscored key).
    # GCA_1: organism name from its report. GCA_2: no report -> accession.
    assert labels == ["Felimare picta", "Brassica rapa", "GCA_2"]


def test_full_jcvi_names_override_every_track(tmp_path):
    (tmp_path / "GCA_1_assembly_report.txt").write_text(REPORT_HEADER)
    labels = resolve_layout_labels(
        _keys(),
        assembly_display_name="Felimare picta",
        custom_names="My Assembly,Species One,Species Two",
        reference_dir=str(tmp_path),
    )
    assert labels == ["My Assembly", "Species One", "Species Two"]


def test_single_jcvi_name_targets_assembly_only(tmp_path):
    (tmp_path / "GCA_1_assembly_report.txt").write_text(REPORT_HEADER)
    labels = resolve_layout_labels(
        _keys(),
        assembly_display_name="Felimare picta",
        custom_names="Just The Assembly",
        reference_dir=str(tmp_path),
    )
    # Single name overrides the assembly; references fall through to report/accession.
    assert labels == ["Just The Assembly", "Brassica rapa", "GCA_2"]


def test_mismatched_count_uses_first_name_for_assembly_only(tmp_path):
    # Count != 1 + N references: matches grid (assembly takes names_list[0],
    # references fall through).
    labels = resolve_layout_labels(
        _keys(),
        assembly_display_name="Felimare picta",
        custom_names="A,B",
        reference_dir=str(tmp_path),
    )
    assert labels == ["A", "GCA_1", "GCA_2"]


def test_assembly_falls_back_to_key_when_no_display_name(tmp_path):
    # When the spaced display name is empty, the caller passes the underscored key.
    labels = resolve_layout_labels(
        ["Felimare_picta"],
        assembly_display_name="Felimare_picta",
        custom_names="",
        reference_dir=str(tmp_path),
    )
    assert labels == ["Felimare_picta"]
