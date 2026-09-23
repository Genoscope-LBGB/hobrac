"""RI threshold CLI propagation and output integration."""

import csv
import sys
from pathlib import Path

import pytest

from hobrac.command_line import get_args
from hobrac.jcvi_synteny import output
from hobrac.main import generate_snakemake_command


@pytest.mark.parametrize("threshold", [None, 1, 10])
def test_hobrac_cli_forwards_ri_threshold(monkeypatch, threshold):
    argv = ["hobrac", "-a", "assembly.fa", "-n", "Species", "-t", "1"]
    if threshold is not None:
        argv += ["--min-ri-alg-genes", str(threshold)]
    monkeypatch.setattr(sys, "argv", argv)
    args = get_args()
    expected = 5 if threshold is None else threshold
    assert args.min_ri_alg_genes == expected
    args.busco_reference_override_paths = {}
    assert f"min_ri_alg_genes={expected} " in generate_snakemake_command(args)


JCVI_ARGS = [
    "jcvi_synteny",
    "--busco_assembly",
    "full_table.tsv",
    "--assembly-fasta",
    "assembly.fa",
    "--busco_references",
    "reference",
    "--accession_order",
    "order.txt",
    "--output_dir",
    "synteny_plots",
]


@pytest.mark.parametrize("threshold", [None, 1, 10])
def test_jcvi_cli_forwards_ri_threshold(monkeypatch, threshold):
    argv = list(JCVI_ARGS)
    if threshold is not None:
        argv += ["--min-ri-alg-genes", str(threshold)]
    monkeypatch.setattr(sys, "argv", argv)
    captured = {}
    monkeypatch.setattr(output, "run", lambda **kwargs: captured.update(kwargs))
    output.main()
    assert captured["min_ri_alg_genes"] == (5 if threshold is None else threshold)


@pytest.mark.parametrize("value", ["0", "-1", "1.5"])
@pytest.mark.parametrize("entry", ["hobrac", "jcvi"])
def test_clis_reject_invalid_threshold(monkeypatch, capsys, value, entry):
    argv = (
        ["hobrac", "-a", "assembly.fa", "-n", "Species", "-t", "1"]
        if entry == "hobrac"
        else list(JCVI_ARGS)
    )
    monkeypatch.setattr(sys, "argv", argv + ["--min-ri-alg-genes", value])
    with pytest.raises(SystemExit) as exc:
        (get_args if entry == "hobrac" else output.main)()
    assert exc.value.code == 2
    assert "--min-ri-alg-genes" in capsys.readouterr().err


@pytest.fixture
def ri_inputs(tmp_path):
    fasta = tmp_path / "assembly.fa"
    fasta.write_text(">chr1\n" + "A" * 2000 + "\n>chr2\n" + "A" * 2000 + "\n")
    busco = tmp_path / "full_table.tsv"
    busco.write_text(
        "".join(f"a{i}\tComplete\tchr1\t{i * 100}\t{i * 100 + 99}\n" for i in range(5))
        + "b0\tComplete\tchr1\t600\t699\n"
        + "b1\tComplete\tchr2\t0\t99\n"
        + "u0\tComplete\tchr1\t700\t799\n"
    )
    colors = tmp_path / "colors.tsv"
    colors.write_text(
        "".join(f"a{i}\t1,2,3\tA\n" for i in range(5))
        + "b0\t4,5,6\tB\nb1\t4,5,6\tB\nmissing\t7,8,9\tC\n"
    )
    order = tmp_path / "order.txt"
    order.write_text("")
    return dict(
        assembly_busco=str(busco),
        assembly_fasta=str(fasta),
        busco_refs=[],
        accession_order=str(order),
        manual_refs="",
        output_dir=str(tmp_path / "out"),
        custom_color_file=str(colors),
        skip_alg=True,
        use_gravity_ordering=False,
        min_busco_genes=2,
    )


def _rows(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.mark.parametrize("threshold,expected_count", [(5, 1), (1, 2), (6, 0)])
def test_run_writes_filtered_summary_and_warnings(ri_inputs, threshold, expected_count):
    result = output.run(**ri_inputs, min_ri_alg_genes=threshold)
    summary = _rows(result["rearrangement_index"])[0]
    assert int(summary["n_algs"]) == expected_count
    detail = _rows(result["rearrangement_index_by_alg"])
    # chr2 fails the sequence filter; the unannotated gene is never counted.
    assert [(row["alg"], row["alg_gene_count"]) for row in detail] == [
        ("A", "5"),
        ("B", "1"),
    ]
    assert detail[0]["chromosome_gene_count"] == "6"
    warnings_path = result["rearrangement_index_warnings"]
    assert Path(warnings_path).name == "rearrangement_index_warnings.tsv"
    warnings = _rows(warnings_path)
    expected_algs = {5: ["B", "C"], 1: ["C"], 6: ["A", "B", "C"]}[threshold]
    assert [row["alg"] for row in warnings] == expected_algs
    assert all(row["min_ri_alg_genes"] == str(threshold) for row in warnings)
    if threshold == 5:
        assert float(summary["rearrangement_index"]) == pytest.approx(1 / 6, abs=1e-6)
    if threshold == 6:
        assert summary["rearrangement_index"] == "nan"


def test_default_run_threshold_is_five(ri_inputs):
    result = output.run(**ri_inputs)
    assert _rows(result["rearrangement_index"])[0]["n_algs"] == "1"
    assert {
        row["min_ri_alg_genes"] for row in _rows(result["rearrangement_index_warnings"])
    } == {"5"}


def test_without_annotations_no_ri_is_produced(ri_inputs):
    result = output.run(**{**ri_inputs, "custom_color_file": ""})
    assert result["rearrangement_index"] == ""
    assert result["rearrangement_index_warnings"] == ""


def test_run_rejects_invalid_threshold_before_reading_inputs(ri_inputs):
    with pytest.raises(ValueError, match="must be at least 1"):
        output.run(**ri_inputs, min_ri_alg_genes=0)
