"""Tests for best-effort chromosome renaming of manual references."""

import gzip
import shlex
import subprocess
import sys
from pathlib import Path
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


@pytest.mark.parametrize(
    "token", ["1", "12", "X", "Y", "Z", "W", "U", "2L", "MT", "Un"]
)
def test_find_chr_name_matches_complete_tokens(token):
    assert find_chr_name(f"AC1 chr{token}, description") == f"chr{token}"
    assert find_chr_name(f"AC1 chromosome: {token}, description") == f"chr{token}"
    assert find_chr_name(f"chr{token}") == f"chr{token}"
    assert find_chr_name(f"chromosome {token}") == f"chr{token}"


@pytest.mark.parametrize("genus", ["Chrysomela", "Chrysochromulina", "Chrysanthemum"])
def test_find_chr_name_does_not_match_genus(genus):
    assert find_chr_name(f"AC1 {genus} unplaced scaffold") is None
    assert find_chr_name(f"AC1 {genus} chromosome 1, complete sequence") == "chr1"
    assert find_chr_name(f"AC1 {genus} chrX, complete sequence") == "chrX"


@pytest.mark.parametrize(
    "header",
    [
        "AC1 chromosome unknown",
        "AC1 chromosome Xylophone",
        "AC1 chromosome 2Long",
        "AC1 microchromosome 1",
        "AC1 _chromosome 1",
        "AC1 mychr1",
        "AC1 _chr1",
        "AC1 chrUn_random",
        "AC1 chr1_random",
        "AC1 chr2Long",
        "AC1 chromosome 1_random",
    ],
)
def test_find_chr_name_rejects_partial_tokens(header):
    assert find_chr_name(header) is None


def test_find_chr_name_matches_case_insensitively():
    assert find_chr_name("AC1 CHRx, description") == "CHRx"
    assert find_chr_name("AC1 CHROMOSOME: mt, description") == "chrmt"


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


def test_rename_reference_ignores_genus_and_preserves_unrecognized_ids(tmp_path):
    src = _write(
        tmp_path / "ref.fa",
        ">AC1 Chrysomela populi chromosome 1, complete sequence\nACGT\n"
        ">AC2 Chrysomela populi unplaced scaffold\nTTTT\n"
        ">chrUn_random unplaced scaffold\nCCCC\n",
    )
    dest = tmp_path / "ref.fna"
    mapping = tmp_path / "ref.chr_rename.tsv"

    result = rename_reference(src, str(dest), str(mapping))

    assert result == [
        ("AC1", "chr1"),
        ("AC2", "AC2"),
        ("chrUn_random", "chrUn_random"),
    ]
    assert dest.read_text() == (
        ">chr1\nACGT\n"
        ">AC2 Chrysomela populi unplaced scaffold\nTTTT\n"
        ">chrUn_random unplaced scaffold\nCCCC\n"
    )
    assert mapping.read_text().splitlines() == [
        "old_name\tnew_name",
        "AC1\tchr1",
        "AC2\tAC2",
        "chrUn_random\tchrUn_random",
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
            result_name = (
                "busco_assembly" if side == "assembly" else "busco_reference_reference"
            )
            tables[side] = (
                output / "busco" / result_name / "run_lineage" / "full_table.tsv"
            )
        else:
            tables[side] = table

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(hobrac_main, "check_dependencies", lambda **kwargs: None)
    monkeypatch.setattr(
        hobrac_main,
        "subprocess",
        SimpleNamespace(
            Popen=lambda *args, **kwargs: SimpleNamespace(returncode=0, wait=lambda: 0)
        ),
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


@pytest.mark.parametrize("reference_count,busco_count", [(0, 1), (1, 2)])
def test_reference_busco_rejects_unmatched_inputs_before_preparation(
    tmp_path, monkeypatch, reference_count, busco_count
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
        "-o",
        str(output),
    ]
    for index in range(reference_count):
        argv.extend(["-r", str(tmp_path / f"reference_{index}.fa")])
    for index in range(busco_count):
        argv.extend(["--busco-reference", str(tmp_path / f"busco_{index}")])
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", argv)

    # Reject the combination before opening inputs or creating output files.
    with pytest.raises(SystemExit) as exit_info:
        hobrac_main.main()
    assert exit_info.value.code == 2
    assert not output.exists()


@pytest.mark.parametrize("reuse_count", [2, 3])
def test_reference_busco_reuses_ordered_prefix(tmp_path, monkeypatch, reuse_count):
    output = tmp_path / "output"
    assembly = tmp_path / "assembly.fa"
    assembly.write_text(">query\n" + "ACGT" * 100 + "\n")
    assembly_table = tmp_path / "assembly_results" / "run_lineage" / "full_table.tsv"
    assembly_table.parent.mkdir(parents=True)
    assembly_table.write_text("123at4751\tComplete\tquery\t10\t100\n")
    argv = [
        "hobrac",
        "-a",
        str(assembly),
        "-n",
        "Test species",
        "-t",
        "1",
        "-e",
        "local",
        "-o",
        str(output),
        "--busco-assembly",
        str(assembly_table),
    ]
    references = []
    # Neither reference names nor BUSCO directory names are alphabetically
    # ordered or paired by name: only the order of the two option lists matters.
    for index, (name, results_name) in enumerate(
        [("zeta", "bundle_two"), ("alpha", "bundle_three"), ("middle", "bundle_one")]
    ):
        fasta = tmp_path / f"{name}.fa.gz"
        with gzip.open(fasta, "wt") as handle:
            handle.write(f">{name}_original chr{index + 1}\n" + "ACGT" * 100 + "\n")
        argv.extend(["-r", str(fasta)])
        sequence_id = f"{name}_original" if index < reuse_count else f"chr{index + 1}"
        table = tmp_path / results_name / "run_lineage" / "full_table.tsv"
        table.parent.mkdir(parents=True)
        table.write_text(f"123at4751\tComplete\t{sequence_id}\t10\t100\n")
        references.append((name, table, sequence_id))

    for _, table, _ in references[:reuse_count]:
        argv.extend(["--busco-reference", str(table)])
    workflow_commands = []

    def defer_workflow(command, **kwargs):
        workflow_commands.append(shlex.split(command))
        return SimpleNamespace(returncode=0, wait=lambda: 0)

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(hobrac_main, "check_dependencies", lambda **kwargs: None)
    monkeypatch.setattr(
        hobrac_main,
        "subprocess",
        SimpleNamespace(Popen=defer_workflow),
    )
    with pytest.raises(SystemExit) as exit_info:
        hobrac_main.main()
    assert exit_info.value.code == 0

    for index, (name, table, sequence_id) in enumerate(references):
        if index < reuse_count:
            table = (
                output
                / "busco"
                / f"busco_reference_{name}"
                / "run_lineage"
                / "full_table.tsv"
            )
        alignment = output / f"alignment_{name}"
        busco_to_paf(
            assembly_table,
            table,
            output / "assembly" / "assembly.fna",
            output / "reference" / f"{name}.fna",
            alignment,
        )
        fields = (alignment / "aln_busco.paf").read_text().strip().split("\t")
        assert [fields[0], fields[5]] == ["query", sequence_id]

    # Inspect the real Snakemake plan for both PAF and JCVI consumers. Supplied
    # result directories must remain inputs, not become targets to recompute.
    (output / "mash").mkdir()
    (output / "mash" / "selected_accessions.txt").write_text("zeta\nalpha\nmiddle\n")
    (output / "busco" / "lineage.txt").write_text("Eukaryota\n")
    (output / "busco" / "datasets.txt").write_text("- eukaryota_odb12\n")
    (output / "busco" / "chosen_dataset.txt").write_text("eukaryota\teukaryota_odb12\n")
    (output / "busco" / "busco_downloads").mkdir()
    command = workflow_commands[0]
    config = command[command.index("--config") + 1 :]
    assert f"assembly_input_fasta={assembly}" in config
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--snakefile",
            hobrac_main.snakefile_path,
            "--dry-run",
            "--printshellcmds",
            "--cores",
            "1",
            "aln/busco_zeta",
            "aln/busco_alpha",
            "aln/busco_middle",
            "synteny_plots/name_mapping.tsv",
            "--config",
            *config,
        ],
        cwd=output,
        capture_output=True,
        text=True,
        timeout=60,
    )
    plan = result.stdout + result.stderr
    assert result.returncode == 0, plan
    assert "--assembly-input-fasta " + shlex.quote(str(assembly)) in plan
    scheduled_fastas = []
    for line in plan.splitlines():
        if line.strip().startswith("busco "):
            tokens = shlex.split(line)
            if "-m" in tokens and tokens[tokens.index("-m") + 1] == "geno":
                scheduled_fastas.append(Path(tokens[tokens.index("-i") + 1]).name)
    assert scheduled_fastas == (["middle.fna"] if reuse_count == 2 else [])
