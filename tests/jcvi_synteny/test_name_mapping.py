"""Name-to-FASTA provenance must follow the exact karyotype track order."""

import csv
from pathlib import Path

import pytest

from hobrac.jcvi_synteny.output import run


@pytest.fixture
def inputs(tmp_path):
    assembly = tmp_path / "prepared.fna"
    assembly.write_text(">chr1\n" + "A" * 200 + "\n")
    table_text = "123at1\tComplete\tchr1\t10\t100\n"
    assembly_busco = tmp_path / "assembly_full_table.tsv"
    assembly_busco.write_text(table_text)
    busco_refs = []
    for key in ["alpha", "zeta"]:
        directory = tmp_path / f"busco_reference_{key}"
        table = directory / "run_lineage" / "full_table.tsv"
        table.parent.mkdir(parents=True)
        table.write_text(table_text)
        busco_refs.append(str(directory))
    order = tmp_path / "selected_accessions.txt"
    order.write_text("zeta\nalpha\n")
    reference_dir = tmp_path / "reference"
    reference_dir.mkdir()
    (reference_dir / "zeta_assembly_report.txt").write_text(
        "# Organism name: Species zeta (common name)\n"
    )
    return dict(
        assembly_busco=str(assembly_busco),
        assembly_fasta=str(assembly),
        busco_refs=busco_refs,
        accession_order=str(order),
        manual_refs="",
        output_dir=str(tmp_path / "synteny_plots"),
        assembly_name="Test_assembly",
        assembly_display_name="Test assembly",
        reference_dir=str(reference_dir),
    )


def read_mapping(result):
    with open(result["name_mapping"], newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.mark.parametrize(
    "names,expected",
    [
        ("", ["Test assembly", "zeta", "alpha"]),
        ("Assembly only", ["Assembly only", "zeta", "alpha"]),
        (
            "My assembly,Second species,First species",
            ["My assembly", "Second species", "First species"],
        ),
    ],
)
def test_manual_mapping_preserves_input_order_and_original_paths(
    inputs, tmp_path, names, expected
):
    # Neither alphabetical BUSCO-directory order nor accession order may
    # override -r ordering. No FASTA read is needed for provenance-only paths.
    Path(inputs["accession_order"]).write_text("alpha\nzeta\n")
    original = tmp_path / "original assembly.fa.gz"
    zeta = tmp_path / "input genomes" / "zeta.fasta.gz"
    alpha = tmp_path / "alpha.fa"
    # Manual references have no assembly reports.
    (Path(inputs["reference_dir"]) / "zeta_assembly_report.txt").unlink()
    result = run(
        **{**inputs, "manual_refs": f"{zeta};{alpha}"},
        assembly_input_fasta=str(original),
        custom_names=names,
    )
    rows = read_mapping(result)
    assert rows == [
        dict(track=str(i), role=role, name=name, identifier=key, fasta=str(fasta))
        for i, role, name, key, fasta in zip(
            [1, 2, 3],
            ["assembly", "reference", "reference"],
            expected,
            ["Test_assembly", "zeta", "alpha"],
            [original, zeta, alpha],
            strict=True,
        )
    ]
    layouts = Path(result["layouts"]).read_text()
    assert "# track_labels\t" + "\t".join(expected) in layouts


@pytest.mark.parametrize("names", ["", "Assembly only", "A,Z,B"])
def test_automatic_mapping_uses_mash_order_and_downloaded_fastas(inputs, names):
    result = run(**inputs, custom_names=names)
    rows = read_mapping(result)
    expected_names = {
        "": ["Test assembly", "Species zeta", "alpha"],
        "Assembly only": ["Assembly only", "Species zeta", "alpha"],
        "A,Z,B": ["A", "Z", "B"],
    }[names]
    assert [row["name"] for row in rows] == expected_names
    assert [row["identifier"] for row in rows] == ["Test_assembly", "zeta", "alpha"]
    assert [row["fasta"] for row in rows] == [
        inputs["assembly_fasta"],
        str(Path(inputs["reference_dir"]) / "zeta.fna"),
        str(Path(inputs["reference_dir"]) / "alpha.fna"),
    ]


def test_mapping_tracks_filtered_species_and_survives_empty_tracks(inputs, tmp_path):
    zeta = tmp_path / "zeta.fa"
    alpha = tmp_path / "alpha.fa"
    missing = tmp_path / "missing.fa"
    result = run(
        **{**inputs, "manual_refs": f"{missing};{zeta};{alpha}"},
        min_busco_genes=30,
    )
    rows = read_mapping(result)
    assert [row["identifier"] for row in rows] == ["Test_assembly", "zeta", "alpha"]
    assert [row["fasta"] for row in rows][1:] == [str(zeta), str(alpha)]
    assert all(not Path(bed).read_text() for bed in result["bed_files"])
