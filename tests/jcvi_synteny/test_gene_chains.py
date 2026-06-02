"""Tests for the per-gene chain table (save_gene_chains)."""

from hobrac.jcvi_synteny.models import BuscoGene
from hobrac.jcvi_synteny.output import save_gene_chains


def _gene(chromosome, start=0, end=100):
    return BuscoGene(busco_id="", chromosome=chromosome, start=start, end=end)


def _read(path):
    with open(path) as f:
        return [line.rstrip("\n").split("\t") for line in f]


def _build(tmp_path):
    # Three species. g1/g2 on chain 0, g3 on chain 1, g4 on no chain.
    # g2 is ABSENT in Bf; g4 is ABSENT in Hl. g1 carries distinct
    # coordinates per species to exercise the *_pos columns.
    species_busco = [
        ("Pm", {"g1": _gene("chr1", 100, 250), "g2": _gene("chr1"),
                "g3": _gene("chr9"), "g4": _gene("chr2")}),
        ("Bf", {"g1": _gene("chr3", 400, 610), "g3": _gene("chr8"),
                "g4": _gene("chr5")}),
        ("Hl", {"g1": _gene("chr20", 700, 820), "g2": _gene("chr20"),
                "g3": _gene("chr7")}),
    ]
    gene_to_chain = {"g1": 0, "g2": 0, "g3": 1, "g4": -1}
    gene_colors = {
        "g1": "#1f77b4",
        "g2": "#1f77b4",
        "g3": "#ff7f0e",
        "g4": "lightgrey",
    }
    out = tmp_path / "gene_chains.tsv"
    save_gene_chains(species_busco, gene_to_chain, gene_colors, str(out))
    return _read(out)


class TestSaveGeneChains:
    def test_header(self, tmp_path):
        rows = _build(tmp_path)
        assert rows[0] == [
            "chain_id", "gene", "Pm", "Bf", "Hl", "color",
            "Pm_pos", "Bf_pos", "Hl_pos",
        ]

    def test_absent_cells(self, tmp_path):
        rows = {r[1]: r for r in _build(tmp_path)[1:]}
        # g2 missing in Bf -> ABSENT; present elsewhere.
        assert rows["g2"][2:5] == ["chr1", "ABSENT", "chr20"]
        # g4 missing in Hl -> ABSENT.
        assert rows["g4"][2:5] == ["chr2", "chr5", "ABSENT"]

    def test_gene_own_chromosome(self, tmp_path):
        rows = {r[1]: r for r in _build(tmp_path)[1:]}
        # Cells are the gene's actual chromosome, which differs per species.
        assert rows["g1"][2:5] == ["chr1", "chr3", "chr20"]

    def test_chain_id_and_color(self, tmp_path):
        rows = {r[1]: r for r in _build(tmp_path)[1:]}
        # color sits right after the chromosome block (index 5 here).
        assert rows["g1"][0] == "0"
        assert rows["g1"][5] == "#1f77b4"
        assert rows["g3"][0] == "1"
        assert rows["g3"][5] == "#ff7f0e"
        # No-chain gene: chain_id -1, lightgrey.
        assert rows["g4"][0] == "-1"
        assert rows["g4"][5] == "lightgrey"

    def test_position_present(self, tmp_path):
        rows = {r[1]: r for r in _build(tmp_path)[1:]}
        # *_pos columns are the three rightmost, per-species start:end.
        assert rows["g1"][6:9] == ["100:250", "400:610", "700:820"]
        # default-coord genes echo their raw ints.
        assert rows["g3"][6:9] == ["0:100", "0:100", "0:100"]

    def test_position_absent(self, tmp_path):
        rows = {r[1]: r for r in _build(tmp_path)[1:]}
        # Missing species -> ABSENT in the *_pos column too.
        assert rows["g2"][6:9] == ["0:100", "ABSENT", "0:100"]
        assert rows["g4"][6:9] == ["0:100", "0:100", "ABSENT"]

    def test_ordering_chain_grouped_grey_last(self, tmp_path):
        order = [r[1] for r in _build(tmp_path)[1:]]
        # Chain 0 (g1,g2) first by id, then chain 1 (g3), grey g4 last.
        assert order == ["g1", "g2", "g3", "g4"]

    def test_empty_writes_header_only(self, tmp_path):
        out = tmp_path / "empty.tsv"
        save_gene_chains([("Pm", {}), ("Bf", {})], {}, {}, str(out))
        rows = _read(out)
        assert rows == [
            ["chain_id", "gene", "Pm", "Bf", "color", "Pm_pos", "Bf_pos"]
        ]

    def test_color_fallback_to_lightgrey(self, tmp_path):
        # Gene missing from gene_colors falls back to lightgrey.
        out = tmp_path / "fb.tsv"
        save_gene_chains([("Pm", {"g9": _gene("chr1", 5, 9)})], {}, {}, str(out))
        rows = {r[1]: r for r in _read(out)[1:]}
        # 1 species: chain_id, gene, Pm, color, Pm_pos
        assert rows["g9"][3] == "lightgrey"
        assert rows["g9"][4] == "5:9"
