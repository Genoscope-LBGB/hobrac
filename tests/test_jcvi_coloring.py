"""Tests for coloring functions and run() branching logic."""

from unittest.mock import patch

import pytest

from hobrac.jcvi_synteny import (
    ALG_PALETTE,
    BuscoGene,
    PairwiseAssociation,
    apply_custom_colors,
    apply_custom_colors_with_algs,
    build_gene_chain_mapping,
    detect_algs_transitive,
    enumerate_chains,
    generate_links_file,
    run,
)

MODULE = "hobrac.jcvi_synteny"


def _gene(chromosome, start=0):
    return BuscoGene(busco_id="", chromosome=chromosome, start=start, end=start + 100)


class TestApplyCustomColors:
    def test_gene_in_custom_gets_custom_color(self):
        sp1 = {"g1": _gene("chr1")}
        sp2 = {"g1": _gene("chrA")}
        result = apply_custom_colors(sp1, sp2, {"g1": "#00ff00"})
        assert result["g1"] == "#00ff00"

    def test_gene_not_in_custom_gets_lightgrey(self):
        sp1 = {"g1": _gene("chr1")}
        sp2 = {"g1": _gene("chrA")}
        result = apply_custom_colors(sp1, sp2, {})
        assert result["g1"] == "lightgrey"

    def test_only_common_genes_included(self):
        sp1 = {"g1": _gene("chr1"), "g2": _gene("chr2")}
        sp2 = {"g1": _gene("chrA"), "g3": _gene("chrC")}
        result = apply_custom_colors(sp1, sp2, {"g1": "#00ff00", "g2": "#0000ff"})
        assert "g1" in result
        assert "g2" not in result
        assert "g3" not in result

    def test_all_genes_in_custom(self):
        sp1 = {"g1": _gene("chr1"), "g2": _gene("chr2")}
        sp2 = {"g1": _gene("chrA"), "g2": _gene("chrB")}
        result = apply_custom_colors(sp1, sp2, {"g1": "#00ff00", "g2": "#0000ff"})
        assert result["g1"] == "#00ff00"
        assert result["g2"] == "#0000ff"

    def test_empty_custom_colors(self):
        sp1 = {"g1": _gene("chr1"), "g2": _gene("chr2")}
        sp2 = {"g1": _gene("chrA"), "g2": _gene("chrB")}
        result = apply_custom_colors(sp1, sp2, {})
        assert result["g1"] == "lightgrey"
        assert result["g2"] == "lightgrey"


class TestApplyCustomColorsWithAlgsEdgeCases:
    @pytest.fixture
    def chain_fixture(self):
        sp1 = {"g1": _gene("chr1"), "g2": _gene("chr3")}
        sp2 = {"g1": _gene("chrA"), "g2": _gene("chrC")}
        gene_to_chain = {"g1": 0, "g2": -1}
        chain_colors = {0: "#ff0000"}
        return sp1, sp2, gene_to_chain, chain_colors

    def test_empty_custom_colors_gives_chain_colors(self, chain_fixture):
        sp1, sp2, gene_to_chain, chain_colors = chain_fixture
        result = apply_custom_colors_with_algs(
            sp1, sp2, gene_to_chain, chain_colors, {}
        )
        assert result["g1"] == chain_colors[0]

    def test_all_genes_in_custom(self, chain_fixture):
        sp1, sp2, gene_to_chain, chain_colors = chain_fixture
        result = apply_custom_colors_with_algs(
            sp1,
            sp2,
            gene_to_chain,
            chain_colors,
            {"g1": "#00ff00", "g2": "#0000ff"},
        )
        assert result["g1"] == "#00ff00"
        assert result["g2"] == "lightgrey"

    def test_no_genes_in_custom(self, chain_fixture):
        sp1, sp2, gene_to_chain, chain_colors = chain_fixture
        result = apply_custom_colors_with_algs(
            sp1, sp2, gene_to_chain, chain_colors, {}
        )
        assert result["g1"] == "#ff0000"
        assert result["g2"] == "lightgrey"


class TestGenerateLinksHideNonSignificant:
    @pytest.fixture
    def links_data(self):
        sp1 = {"g1": _gene("chr1", 1000), "g2": _gene("chr3", 1000)}
        sp2 = {"g1": _gene("chrA", 1000), "g2": _gene("chrC", 1000)}
        gene_colors = {"g1": "#ff0000", "g2": "lightgrey"}
        return sp1, sp2, gene_colors

    def test_lightgrey_excluded_colored_included(self, tmp_path, links_data):
        sp1, sp2, gene_colors = links_data
        path = str(tmp_path / "links.simple")
        generate_links_file(
            sp1, sp2, gene_colors, "sp1", "sp2", path, hide_non_significant=True
        )
        with open(path) as f:
            content = f.read()
        assert "#ff0000" in content
        assert "lightgrey" not in content

    def test_all_blocks_when_hide_false(self, tmp_path, links_data):
        sp1, sp2, gene_colors = links_data
        path = str(tmp_path / "links.simple")
        generate_links_file(
            sp1, sp2, gene_colors, "sp1", "sp2", path, hide_non_significant=False
        )
        with open(path) as f:
            content = f.read()
        assert "#ff0000" in content
        assert "lightgrey" in content


class TestRunBranching:
    @pytest.fixture
    def run_mocks(self):
        busco_data = {
            "g1": BuscoGene(busco_id="g1", chromosome="chr1", start=0, end=100),
        }
        mock_specs = {
            "read_fasta_sizes": {"return_value": {}},
            "read_busco_tsv": {"return_value": busco_data},
            "glob.glob": {"return_value": ["/fake/full_table.tsv"]},
            "generate_bed_file": {},
            "generate_links_file": {},
            "generate_seqids_file": {},
            "generate_layouts_file": {},
            "save_alg_associations": {},
            "detect_algs_transitive": {"return_value": ([], [], {}, {})},
            "build_alg_association_list": {"return_value": []},
            "apply_custom_colors": {"return_value": {"g1": "lightgrey"}},
            "apply_custom_colors_with_algs": {"return_value": {"g1": "lightgrey"}},
            "parse_custom_colors": {"return_value": {"g1": "#00ff00"}},
        }
        patchers = []
        mocks = {}
        for name, kwargs in mock_specs.items():
            p = patch(f"{MODULE}.{name}", **kwargs)
            mocks[name] = p.start()
            patchers.append(p)
        yield mocks
        for p in patchers:
            p.stop()

    def _call_run(self, tmp_path, custom_color_file="", skip_alg=False):
        run(
            assembly_busco="/fake/assembly_busco",
            assembly_fasta="/fake/assembly.fasta",
            busco_refs=["/fake/busco_reference_ref1"],
            accession_order="",
            manual_refs="/fake/ref1.fasta",
            output_dir=str(tmp_path / "output"),
            custom_color_file=custom_color_file,
            skip_alg=skip_alg,
        )

    def test_custom_with_skip_alg(self, tmp_path, run_mocks):
        self._call_run(tmp_path, custom_color_file="/fake/colors.tsv", skip_alg=True)

        run_mocks["apply_custom_colors"].assert_called()
        run_mocks["detect_algs_transitive"].assert_not_called()
        run_mocks["apply_custom_colors_with_algs"].assert_not_called()

    def test_custom_without_skip_alg(self, tmp_path, run_mocks):
        self._call_run(tmp_path, custom_color_file="/fake/colors.tsv", skip_alg=False)

        run_mocks["apply_custom_colors_with_algs"].assert_called()
        run_mocks["detect_algs_transitive"].assert_called()
        run_mocks["apply_custom_colors"].assert_not_called()

    def test_no_custom_colors(self, tmp_path, run_mocks):
        self._call_run(tmp_path, custom_color_file="", skip_alg=False)

        run_mocks["apply_custom_colors_with_algs"].assert_called()
        run_mocks["detect_algs_transitive"].assert_called()
        run_mocks["apply_custom_colors"].assert_not_called()


def _assoc(sp1, sp2, chr1, chr2):
    return PairwiseAssociation(sp1, sp2, chr1, chr2, p_value=0.001, gene_count=10)


class TestEnumerateChains:
    def test_empty_input(self):
        assert enumerate_chains([]) == []

    def test_two_species_single_pair(self):
        assocs = [_assoc("A", "B", "A1", "B1")]
        chains = enumerate_chains(assocs)
        assert chains == [[("A", "A1"), ("B", "B1")]]

    def test_two_species_multiple_pairs(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("A", "B", "A2", "B2"),
        ]
        chains = enumerate_chains(assocs)
        assert len(chains) == 2
        assert [("A", "A1"), ("B", "B1")] in chains
        assert [("A", "A2"), ("B", "B2")] in chains

    def test_three_species_linear(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("B", "C", "B1", "C1"),
        ]
        chains = enumerate_chains(assocs)
        assert chains == [[("A", "A1"), ("B", "B1"), ("C", "C1")]]

    def test_three_species_downstream_branching(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("B", "C", "B1", "C2"),
            _assoc("B", "C", "B1", "C3"),
        ]
        chains = enumerate_chains(assocs)
        assert len(chains) == 2
        assert [("A", "A1"), ("B", "B1"), ("C", "C2")] in chains
        assert [("A", "A1"), ("B", "B1"), ("C", "C3")] in chains

    def test_three_species_upstream_merging(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("A", "B", "A2", "B1"),
            _assoc("B", "C", "B1", "C2"),
        ]
        chains = enumerate_chains(assocs)
        assert len(chains) == 2
        assert [("A", "A1"), ("B", "B1"), ("C", "C2")] in chains
        assert [("A", "A2"), ("B", "B1"), ("C", "C2")] in chains

    def test_three_species_dead_end(self):
        assocs = [_assoc("A", "B", "A1", "B1")]
        chains = enumerate_chains(assocs)
        assert chains == [[("A", "A1"), ("B", "B1")]]

    def test_orphan_chain_later_species(self):
        assocs = [_assoc("B", "C", "B3", "C4")]
        chains = enumerate_chains(assocs)
        assert chains == [[("B", "B3"), ("C", "C4")]]

    def test_deterministic_ordering(self):
        assocs = [
            _assoc("A", "B", "A2", "B2"),
            _assoc("A", "B", "A1", "B1"),
        ]
        chains = enumerate_chains(assocs)
        assert chains[0] == [("A", "A1"), ("B", "B1")]
        assert chains[1] == [("A", "A2"), ("B", "B2")]

    def test_mixed_dead_end_and_orphan(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("B", "C", "B3", "C4"),
        ]
        chains = enumerate_chains(assocs)
        assert len(chains) == 2
        assert [("A", "A1"), ("B", "B1")] in chains
        assert [("B", "B3"), ("C", "C4")] in chains

    def test_four_species_spanning_chain(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("B", "C", "B1", "C2"),
            _assoc("C", "D", "C2", "D3"),
        ]
        chains = enumerate_chains(assocs)
        assert chains == [[("A", "A1"), ("B", "B1"), ("C", "C2"), ("D", "D3")]]

    def test_diamond_converging_at_last_species(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("A", "B", "A1", "B2"),
            _assoc("B", "C", "B1", "C2"),
            _assoc("B", "C", "B2", "C2"),
        ]
        chains = enumerate_chains(assocs)
        assert len(chains) == 2
        assert [("A", "A1"), ("B", "B1"), ("C", "C2")] in chains
        assert [("A", "A1"), ("B", "B2"), ("C", "C2")] in chains


class TestBuildGeneChainMapping:
    def test_exact_match_three_species(self):
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B1")}),
            ("C", {"g1": _gene("C2")}),
        ]
        chains = [[("A", "A1"), ("B", "B1"), ("C", "C2")]]
        result = build_gene_chain_mapping(species_busco, chains)
        assert result["g1"] == 0

    def test_no_match_gives_minus_one(self):
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B9")}),
        ]
        chains = [[("A", "A1"), ("B", "B1")]]
        result = build_gene_chain_mapping(species_busco, chains)
        assert result["g1"] == -1

    def test_partial_match_gene_absent_from_one_species(self):
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B1")}),
            ("C", {}),
        ]
        chains = [[("A", "A1"), ("B", "B1"), ("C", "C2")]]
        result = build_gene_chain_mapping(species_busco, chains)
        assert result["g1"] == 0

    def test_ambiguous_picks_most_populated_chain(self):
        # g2, g4 exactly match chain 0; g3 exactly matches chain 1
        # g1 absent from C → matches both → picks chain 0 (2 genes vs 1)
        species_busco = [
            ("A", {k: _gene("A1") for k in ("g1", "g2", "g3", "g4")}),
            ("B", {k: _gene("B1") for k in ("g1", "g2", "g3", "g4")}),
            ("C", {"g2": _gene("C2"), "g3": _gene("C3"), "g4": _gene("C2")}),
        ]
        chains = [
            [("A", "A1"), ("B", "B1"), ("C", "C2")],
            [("A", "A1"), ("B", "B1"), ("C", "C3")],
        ]
        result = build_gene_chain_mapping(species_busco, chains)
        assert result["g2"] == 0
        assert result["g4"] == 0
        assert result["g3"] == 1
        assert result["g1"] == 0

    def test_ambiguous_tiebreak_lowest_id(self):
        # Equal population → tie-break by lowest chain ID
        species_busco = [
            ("A", {k: _gene("A1") for k in ("g1", "g2", "g3")}),
            ("B", {k: _gene("B1") for k in ("g1", "g2", "g3")}),
            ("C", {"g2": _gene("C2"), "g3": _gene("C3")}),
        ]
        chains = [
            [("A", "A1"), ("B", "B1"), ("C", "C2")],
            [("A", "A1"), ("B", "B1"), ("C", "C3")],
        ]
        result = build_gene_chain_mapping(species_busco, chains)
        assert result["g1"] == 0

    def test_two_species_baseline(self):
        species_busco = [
            ("A", {"g1": _gene("A1"), "g2": _gene("A2")}),
            ("B", {"g1": _gene("B1"), "g2": _gene("B2")}),
        ]
        chains = [
            [("A", "A1"), ("B", "B1")],
            [("A", "A2"), ("B", "B2")],
        ]
        result = build_gene_chain_mapping(species_busco, chains)
        assert result["g1"] == 0
        assert result["g2"] == 1

    def test_four_species(self):
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B1")}),
            ("C", {"g1": _gene("C2")}),
            ("D", {"g1": _gene("D3")}),
        ]
        chains = [[("A", "A1"), ("B", "B1"), ("C", "C2"), ("D", "D3")]]
        result = build_gene_chain_mapping(species_busco, chains)
        assert result["g1"] == 0

    def test_gene_in_single_species_excluded(self):
        species_busco = [
            ("A", {"g1": _gene("A1"), "lone": _gene("A9")}),
            ("B", {"g1": _gene("B1")}),
        ]
        chains = [[("A", "A1"), ("B", "B1")]]
        result = build_gene_chain_mapping(species_busco, chains)
        assert "lone" not in result

    def test_empty_chains_all_minus_one(self):
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B1")}),
        ]
        result = build_gene_chain_mapping(species_busco, [])
        assert result["g1"] == -1

    def test_empty_species_returns_empty(self):
        result = build_gene_chain_mapping([], [[("A", "A1")]])
        assert result == {}


class TestDetectAlgsTransitive:
    @patch(f"{MODULE}.detect_algs_pairwise_raw")
    def test_two_species_returns_gene_to_chain(self, mock_pairwise):
        sp1 = {"g1": _gene("A1"), "g2": _gene("A2")}
        sp2 = {"g1": _gene("B1"), "g2": _gene("B2")}
        mock_pairwise.return_value = [
            _assoc("sp1", "sp2", "A1", "B1"),
            _assoc("sp1", "sp2", "A2", "B2"),
        ]

        assocs, chains, gene_to_chain, chain_colors = detect_algs_transitive(
            [("sp1", sp1), ("sp2", sp2)]
        )

        assert len(assocs) == 2
        assert len(chains) == 2
        assert isinstance(gene_to_chain, dict)
        assert all(isinstance(v, int) for v in gene_to_chain.values())
        assert gene_to_chain["g1"] != gene_to_chain["g2"]
        assert gene_to_chain["g1"] in chain_colors
        assert gene_to_chain["g2"] in chain_colors

    @patch(f"{MODULE}.detect_algs_pairwise_raw")
    def test_each_chain_gets_palette_color(self, mock_pairwise):
        sp1 = {"g1": _gene("A1"), "g2": _gene("A2")}
        sp2 = {"g1": _gene("B1"), "g2": _gene("B2")}
        mock_pairwise.return_value = [
            _assoc("sp1", "sp2", "A1", "B1"),
            _assoc("sp1", "sp2", "A2", "B2"),
        ]

        _, _, _, chain_colors = detect_algs_transitive([("sp1", sp1), ("sp2", sp2)])

        for chain_id, color in chain_colors.items():
            assert color == ALG_PALETTE[chain_id % len(ALG_PALETTE)]

    @patch(f"{MODULE}.detect_algs_pairwise_raw")
    def test_three_species_branching_distinct_chains(self, mock_pairwise):
        sp1 = {"g1": _gene("A1"), "g2": _gene("A1")}
        sp2 = {"g1": _gene("B1"), "g2": _gene("B1")}
        sp3 = {"g1": _gene("C1"), "g2": _gene("C2")}
        mock_pairwise.side_effect = [
            [_assoc("sp1", "sp2", "A1", "B1")],
            [_assoc("sp2", "sp3", "B1", "C1"), _assoc("sp2", "sp3", "B1", "C2")],
        ]

        _, _, gene_to_chain, chain_colors = detect_algs_transitive(
            [("sp1", sp1), ("sp2", sp2), ("sp3", sp3)]
        )

        assert gene_to_chain["g1"] != gene_to_chain["g2"]
        assert len(chain_colors) == 2

    @patch(f"{MODULE}.detect_algs_pairwise_raw")
    def test_no_associations_returns_empty(self, mock_pairwise):
        sp1 = {"g1": _gene("A1")}
        sp2 = {"g1": _gene("B1")}
        mock_pairwise.return_value = []

        assocs, chains, gene_to_chain, chain_colors = detect_algs_transitive(
            [("sp1", sp1), ("sp2", sp2)]
        )

        assert assocs == []
        assert chains == []
        assert chain_colors == {}
        assert all(v == -1 for v in gene_to_chain.values())

    @patch(f"{MODULE}.detect_algs_pairwise_raw")
    def test_single_species_returns_empty(self, mock_pairwise):
        sp1 = {"g1": _gene("A1")}

        assocs, chains, gene_to_chain, chain_colors = detect_algs_transitive(
            [("sp1", sp1)]
        )

        mock_pairwise.assert_not_called()
        assert assocs == []
        assert chains == []
        assert gene_to_chain == {}
        assert chain_colors == {}


class TestThreeSpeciesIntegration:
    """Integration tests: BUSCO data → detect_algs_transitive → coloring."""

    @pytest.fixture
    def branching_species(self):
        """A1-B1-C2 and A1-B2-C3: genes differ by B chromosome linkage."""
        sp1 = {
            "g_b1": _gene("A1"),
            "g_b2": _gene("A1"),
            "g_nosig": _gene("A9"),
        }
        sp2 = {
            "g_b1": _gene("B1"),
            "g_b2": _gene("B2"),
            "g_nosig": _gene("B9"),
        }
        sp3 = {
            "g_b1": _gene("C2"),
            "g_b2": _gene("C3"),
            "g_nosig": _gene("C9"),
        }
        return [("sp1", sp1), ("sp2", sp2), ("sp3", sp3)]

    @pytest.fixture
    def branching_pairwise(self):
        return [
            [_assoc("sp1", "sp2", "A1", "B1"), _assoc("sp1", "sp2", "A1", "B2")],
            [_assoc("sp2", "sp3", "B1", "C2"), _assoc("sp2", "sp3", "B2", "C3")],
        ]

    @pytest.fixture
    def linear_species(self):
        """A1-B1-C2: single linear chain."""
        sp1 = {
            "g_sig": _gene("A1"),
            "g_nosig": _gene("A9"),
        }
        sp2 = {
            "g_sig": _gene("B1"),
            "g_nosig": _gene("B9"),
        }
        sp3 = {
            "g_sig": _gene("C2"),
            "g_nosig": _gene("C9"),
        }
        return [("sp1", sp1), ("sp2", sp2), ("sp3", sp3)]

    @pytest.fixture
    def linear_pairwise(self):
        return [
            [_assoc("sp1", "sp2", "A1", "B1")],
            [_assoc("sp2", "sp3", "B1", "C2")],
        ]

    @patch(f"{MODULE}.detect_algs_pairwise_raw")
    def test_branching_genes_get_different_colors(
        self, mock_pairwise, branching_species, branching_pairwise
    ):
        mock_pairwise.side_effect = branching_pairwise

        _, _, gene_to_chain, chain_colors = detect_algs_transitive(branching_species)

        assert gene_to_chain["g_b1"] != gene_to_chain["g_b2"]
        assert gene_to_chain["g_b1"] >= 0
        assert gene_to_chain["g_b2"] >= 0
        assert gene_to_chain["g_nosig"] == -1

        sp1_busco = branching_species[0][1]
        sp2_busco = branching_species[1][1]
        colors_ab = apply_custom_colors_with_algs(
            sp1_busco, sp2_busco, gene_to_chain, chain_colors, {}
        )
        assert colors_ab["g_b1"] != colors_ab["g_b2"]
        assert colors_ab["g_b1"] != "lightgrey"
        assert colors_ab["g_b2"] != "lightgrey"
        assert colors_ab["g_nosig"] == "lightgrey"

    @patch(f"{MODULE}.detect_algs_pairwise_raw")
    def test_linear_chain_same_color_across_pairs(
        self, mock_pairwise, linear_species, linear_pairwise
    ):
        mock_pairwise.side_effect = linear_pairwise

        _, _, gene_to_chain, chain_colors = detect_algs_transitive(linear_species)

        colors_ab = apply_custom_colors_with_algs(
            linear_species[0][1], linear_species[1][1], gene_to_chain, chain_colors, {}
        )
        colors_bc = apply_custom_colors_with_algs(
            linear_species[1][1], linear_species[2][1], gene_to_chain, chain_colors, {}
        )

        assert colors_ab["g_sig"] == colors_bc["g_sig"]
        assert colors_ab["g_sig"] != "lightgrey"

    @patch(f"{MODULE}.detect_algs_pairwise_raw")
    def test_linear_chain_with_custom_colors(
        self, mock_pairwise, linear_species, linear_pairwise
    ):
        mock_pairwise.side_effect = linear_pairwise

        _, _, gene_to_chain, chain_colors = detect_algs_transitive(linear_species)

        custom = {"g_sig": "#00ff00"}
        sp1_busco = linear_species[0][1]
        sp2_busco = linear_species[1][1]

        colors = apply_custom_colors_with_algs(
            sp1_busco, sp2_busco, gene_to_chain, chain_colors, custom
        )

        assert colors["g_sig"] == "#00ff00"
        assert colors["g_nosig"] == "lightgrey"


class TestTwoSpeciesBackwardCompat:
    """Verify 2-species chain-based pipeline produces correct coloring."""

    @patch(f"{MODULE}.detect_algs_pairwise_raw")
    def test_two_species_coloring_matches_chain_assignment(self, mock_pairwise):
        sp1 = {
            "g1": _gene("A1"),
            "g2": _gene("A2"),
            "g3": _gene("A1"),
            "g_nosig": _gene("A9"),
        }
        sp2 = {
            "g1": _gene("B1"),
            "g2": _gene("B2"),
            "g3": _gene("B1"),
            "g_nosig": _gene("B9"),
        }
        mock_pairwise.return_value = [
            _assoc("sp1", "sp2", "A1", "B1"),
            _assoc("sp1", "sp2", "A2", "B2"),
        ]

        _, chains, gene_to_chain, chain_colors = detect_algs_transitive(
            [("sp1", sp1), ("sp2", sp2)]
        )

        colors = apply_custom_colors_with_algs(
            sp1, sp2, gene_to_chain, chain_colors, {}
        )

        assert colors["g1"] == colors["g3"]
        assert colors["g1"] != colors["g2"]
        assert colors["g1"] in ALG_PALETTE
        assert colors["g2"] in ALG_PALETTE
        assert colors["g_nosig"] == "lightgrey"
        assert len(chains) == 2
        for cid, color in chain_colors.items():
            assert color == ALG_PALETTE[cid % len(ALG_PALETTE)]
