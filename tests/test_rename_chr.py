"""Tests for best-effort chromosome renaming of manual references."""

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
        ">CM000663.2 Homo sapiens chr1, GRCh38\n"
        "ACGT\n"
        ">scaffold_42 unplaced\n"
        "TTTT\n",
    )
    dest = tmp_path / "ref.fna"
    mapping = tmp_path / "ref.chr_rename.tsv"

    result = rename_reference(src, str(dest), str(mapping))

    assert result == [("CM000663.2", "chr1"), ("scaffold_42", "scaffold_42")]

    # Renamed header is replaced; unchanged header is kept verbatim, and the
    # sequence data is preserved.
    assert dest.read_text() == (
        ">chr1\n"
        "ACGT\n"
        ">scaffold_42 unplaced\n"
        "TTTT\n"
    )

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
