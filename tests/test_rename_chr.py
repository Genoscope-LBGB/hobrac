"""Tests for best-effort chromosome renaming of manual references."""

import gzip
import sys
from types import SimpleNamespace

import pytest

import hobrac.main as hobrac_main
from hobrac.busco_to_paf import run as busco_to_paf
from hobrac.rename_chr import find_chr_name, rename_reference


def test_find_chr_name_matches_common_tokens():
    assert find_chr_name("chr1 some description") == "chr1"
    assert find_chr_name("scaffold chrX more") == "chrX"
    assert find_chr_name("chr2L drosophila arm") == "chr2L"
    assert find_chr_name("chrMT mitochondrion") == "chrMT"


def test_find_chr_name_matches_descriptive_chromosome():
    # GenBank / ENA style descriptions are normalized to chr<token>.
    assert (
        find_chr_name(
            "CM090417.1 Ctenoides ales isolate KM-2024 chromosome 1,"
            " whole genome shotgun sequence"
        )
        == "chr1"
    )
    assert (
        find_chr_name("LR736838.1 Pecten maximus genome assembly, chromosome: 1")
        == "chr1"
    )
    assert (
        find_chr_name("OZ121646.1 Venus verrucosa genome assembly, chromosome: 4")
        == "chr4"
    )
    assert find_chr_name("AC1 genome assembly, chromosome: 2L") == "chr2L"
    assert find_chr_name("AC1 genome assembly, chromosome X") == "chrX"


def test_find_chr_name_ignores_assembly_name_and_plurals():
    # The Ensembl "chromosome:ASSEMBLY:NAME" form must not be misread, and
    # plurals must be ignored.
    assert find_chr_name("1 dna:chromosome chromosome:GRCh38:1") is None
    assert find_chr_name("scaffold spanning 3 chromosomes") is None


def test_find_chr_name_no_match():
    assert find_chr_name("scaffold_123 unplaced") is None


def _write(path, text):
    path.write_text(text)
    return str(path)


def test_rename_reference_renames_and_writes_mapping(tmp_path):
    src = _write(
        tmp_path / "ref.fa",
        ">CM000663.2 Homo sapiens chr1, GRCh38\nACGT\n>scaffold_42 unplaced\nTTTT\n",
    )
    dest = tmp_path / "ref.fna"
    mapping = tmp_path / "ref.chr_rename.tsv"

    result = rename_reference(src, str(dest), str(mapping))

    assert result == [("CM000663.2", "chr1"), ("scaffold_42", "scaffold_42")]

    # Renamed header is replaced; unchanged header is kept verbatim, and the
    # sequence data is preserved.
    assert dest.read_text() == (">chr1\nACGT\n>scaffold_42 unplaced\nTTTT\n")

    assert mapping.read_text().splitlines() == [
        "old_name\tnew_name",
        "CM000663.2\tchr1",
        "scaffold_42\tscaffold_42",
    ]


def test_rename_reference_avoids_collisions(tmp_path):
    src = _write(
        tmp_path / "ref.fa",
        ">a chr1\nAA\n>b chr1\nCC\n",
    )
    dest = tmp_path / "ref.fna"
    mapping = tmp_path / "map.tsv"

    result = rename_reference(src, str(dest), str(mapping))

    # First sequence wins the chr1 name; the second keeps its original id.
    assert result == [("a", "chr1"), ("b", "b")]


@pytest.mark.parametrize(
    "reuse_assembly,reuse_reference",
    [(False, False), (False, True), (True, False), (True, True)],
)
def test_busco_reuse_preserves_matching_ids_per_input(
    tmp_path, monkeypatch, reuse_assembly, reuse_reference
):
    output = tmp_path / "output"
    argv = [
        "hobrac",
        "-n",
        "Test species",
        "-t",
        "1",
        "-e",
        "local",
        "-o",
        str(output),
    ]
    tables = {}
    expected_ids = []
    for side, option, reuse, original_id, chromosome in [
        ("assembly", "-a", reuse_assembly, "scaffoldA", "chr1"),
        ("reference", "-r", reuse_reference, "NC_000002.1", "chr2"),
    ]:
        fasta = tmp_path / f"{side}.fa.gz"
        with gzip.open(fasta, "wt") as handle:
            handle.write(f">{original_id} {chromosome}\n" + "ACGT" * 100 + "\n")
        argv.extend([option, str(fasta)])

        # Reused results name the original FASTA; newly computed results name
        # the prepared FASTA. Both must remain usable by downstream conversion.
        sequence_id = original_id if reuse else chromosome
        expected_ids.append(sequence_id)
        table = tmp_path / f"busco_{side}" / "run_lineage" / "full_table.tsv"
        table.parent.mkdir(parents=True)
        table.write_text(f"123at4751\tComplete\t{sequence_id}\t10\t100\n")
        if reuse:
            argv.extend([f"--busco-{side}", str(table)])
            tables[side] = (
                output / "busco" / f"busco_{side}" / "run_lineage" / "full_table.tsv"
            )
        else:
            tables[side] = table

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(hobrac_main, "check_dependencies", lambda **kwargs: None)
    monkeypatch.setattr(
        hobrac_main.subprocess,
        "Popen",
        lambda *args, **kwargs: SimpleNamespace(returncode=0, wait=lambda: 0),
    )
    with pytest.raises(SystemExit) as exit_info:
        hobrac_main.main()
    assert exit_info.value.code == 0

    alignment = output / "alignment"
    busco_to_paf(
        tables["assembly"],
        tables["reference"],
        output / "assembly" / "assembly.fna",
        output / "reference" / "reference.fna",
        alignment,
    )
    fields = (alignment / "aln_busco.paf").read_text().strip().split("\t")
    assert [fields[0], fields[5]] == expected_ids


@pytest.mark.parametrize("references", [(), ("reference_1.fa", "reference_2.fa")])
def test_busco_reference_requires_one_manual_reference_before_preparation(
    tmp_path, monkeypatch, references
):
    output = tmp_path / "output"
    argv = [
        "hobrac",
        "-a",
        str(tmp_path / "assembly.fa"),
        "-n",
        "Test species",
        "-t",
        "1",
        "--busco-reference",
        str(tmp_path / "busco_reference"),
        "-o",
        str(output),
    ]
    for reference in references:
        argv.extend(["-r", str(tmp_path / reference)])
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", argv)

    # Reject the combination before opening inputs or creating output files.
    with pytest.raises(SystemExit) as exit_info:
        hobrac_main.main()
    assert exit_info.value.code == 2
    assert not output.exists()
