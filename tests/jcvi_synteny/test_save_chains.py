"""Tests for save_chains: algs.tsv wide-format output."""

from hobrac.jcvi_synteny.output import save_chains


def _read(path):
    with open(path) as f:
        return [line.rstrip("\n").split("\t") for line in f]


class TestSaveChains:
    def test_writes_header_and_wide_rows_with_dash_for_uncovered(self, tmp_path):
        chains = [
            [("assembly", "chr1"), ("ref_A", "chrX"), ("ref_B", "chr4")],
            [("ref_A", "chr3"), ("ref_B", "chr2")],
        ]
        chain_colors = {0: "#e6194b", 1: "#3cb44b"}
        gene_to_chain = {
            "g1": 0, "g2": 0, "g3": 0,
            "g4": 1, "g5": 1,
            "g6": -1,
        }
        species_names = ["assembly", "ref_A", "ref_B"]
        out = tmp_path / "algs.tsv"

        save_chains(chains, chain_colors, gene_to_chain, species_names, str(out))

        rows = _read(out)
        assert rows[0] == [
            "chain_id", "color", "n_genes", "assembly", "ref_A", "ref_B",
        ]
        assert rows[1] == ["0", "#e6194b", "3", "chr1", "chrX", "chr4"]
        assert rows[2] == ["1", "#3cb44b", "2", "-", "chr3", "chr2"]

    def test_header_only_when_no_chains(self, tmp_path):
        out = tmp_path / "algs.tsv"
        save_chains([], {}, {}, ["assembly", "ref_A"], str(out))

        rows = _read(out)
        assert len(rows) == 1
        assert rows[0] == ["chain_id", "color", "n_genes", "assembly", "ref_A"]

    def test_unassigned_genes_do_not_inflate_counts(self, tmp_path):
        chains = [[("assembly", "chr1"), ("ref_A", "chrX")]]
        chain_colors = {0: "#e6194b"}
        gene_to_chain = {"g1": 0, "g2": -1, "g3": -1}
        out = tmp_path / "algs.tsv"

        save_chains(
            chains, chain_colors, gene_to_chain, ["assembly", "ref_A"], str(out),
        )

        rows = _read(out)
        assert rows[1][2] == "1"

    def test_row_order_matches_chain_id(self, tmp_path):
        chains = [
            [("a", "c0"), ("b", "c0")],
            [("a", "c1"), ("b", "c1")],
            [("a", "c2"), ("b", "c2")],
        ]
        chain_colors = {0: "#111111", 1: "#222222", 2: "#333333"}
        out = tmp_path / "algs.tsv"

        save_chains(chains, chain_colors, {}, ["a", "b"], str(out))

        rows = _read(out)
        assert [r[0] for r in rows[1:]] == ["0", "1", "2"]
