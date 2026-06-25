"""Tests for coloring functions and run() branching logic."""

from unittest.mock import patch

import pytest

from hobrac.jcvi_synteny.chains import build_gene_chain_mapping, enumerate_chains
from hobrac.jcvi_synteny.coloring import (
    apply_custom_colors,
    apply_custom_colors_with_algs,
)
from hobrac.jcvi_synteny.models import ALG_PALETTE, BuscoGene, PairwiseAssociation
from hobrac.jcvi_synteny.output import generate_links_file, run
from hobrac.jcvi_synteny.statistics import detect_algs_transitive

MODULE = "hobrac.jcvi_synteny.output"
STATS_MODULE = "hobrac.jcvi_synteny.statistics"


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
        chains = [[("sp1", "chr1"), ("sp2", "chrA")]]
        return sp1, sp2, gene_to_chain, chain_colors, chains

    def test_empty_custom_colors_gives_chain_colors(self, chain_fixture):
        sp1, sp2, gene_to_chain, chain_colors, chains = chain_fixture
        result = apply_custom_colors_with_algs(
            sp1, sp2, gene_to_chain, chain_colors, {},
            chains=chains, sp1_name="sp1", sp2_name="sp2",
        )
        assert result["g1"] == chain_colors[0]

    def test_all_genes_in_custom(self, chain_fixture):
        sp1, sp2, gene_to_chain, chain_colors, chains = chain_fixture
        result = apply_custom_colors_with_algs(
            sp1,
            sp2,
            gene_to_chain,
            chain_colors,
            {"g1": "#00ff00", "g2": "#0000ff"},
            chains=chains, sp1_name="sp1", sp2_name="sp2",
        )
        # Custom colors are gated by chain membership: g1 is on a chain and
        # listed → custom; g2 is off-chain → lightgrey even though listed.
        assert result["g1"] == "#00ff00"
        assert result["g2"] == "lightgrey"

    def test_no_genes_in_custom(self, chain_fixture):
        sp1, sp2, gene_to_chain, chain_colors, chains = chain_fixture
        result = apply_custom_colors_with_algs(
            sp1, sp2, gene_to_chain, chain_colors, {},
            chains=chains, sp1_name="sp1", sp2_name="sp2",
        )
        assert result["g1"] == "#ff0000"
        assert result["g2"] == "lightgrey"

    def test_gene_grey_when_chain_does_not_cover_pair(self):
        """A gene on chain [(B,B1),(C,C1)] must be lightgrey in the A-B pair
        because the chain does not cover the A→B edge."""
        sp_a = {"g1": _gene("A1")}
        sp_b = {"g1": _gene("B1")}
        sp_c = {"g1": _gene("C1")}
        # Chain only covers B→C, not A→B
        chains = [[("B", "B1"), ("C", "C1")]]
        gene_to_chain = {"g1": 0}
        chain_colors = {0: "#ff0000"}

        # A-B pair: chain doesn't cover this edge → lightgrey
        colors_ab = apply_custom_colors_with_algs(
            sp_a, sp_b, gene_to_chain, chain_colors, {},
            chains=chains, sp1_name="A", sp2_name="B",
        )
        assert colors_ab["g1"] == "lightgrey"

        # B-C pair: chain covers this edge → colored
        colors_bc = apply_custom_colors_with_algs(
            sp_b, sp_c, gene_to_chain, chain_colors, {},
            chains=chains, sp1_name="B", sp2_name="C",
        )
        assert colors_bc["g1"] == "#ff0000"

    def test_chain_covers_pair_canonical_order(self):
        """Chain in canonical (lex) order must still cover a pair whose
        species_busco order is non-lex."""
        sp_z = {"g1": _gene("Z1")}
        sp_a = {"g1": _gene("A1")}
        # Chain in canonical order: A before Z
        chains = [[("A", "A1"), ("Z", "Z1")]]
        gene_to_chain = {"g1": 0}
        chain_colors = {0: "#ff0000"}

        # species_busco order is Z, A (non-lex)
        colors = apply_custom_colors_with_algs(
            sp_z, sp_a, gene_to_chain, chain_colors, {},
            chains=chains, sp1_name="Z", sp2_name="A",
        )
        assert colors["g1"] == "#ff0000"


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


def _parse_simple_file(path):
    """Return a dict {busco_id: color} parsed from a .simple file."""
    result = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            assert len(parts) == 6, f"expected 6 columns, got {len(parts)}: {parts!r}"
            first = parts[0]
            assert "*" in first, f"missing color prefix in {first!r}"
            color, gene1_start = first.split("*", 1)
            # Format: "{sp}_{busco_id}", and start == end for single-gene lines
            assert (
                parts[1] == gene1_start
            ), "start1 != end1 (blocks must be single-gene)"
            assert parts[2] == parts[3], "start2 != end2 (blocks must be single-gene)"
            assert parts[4] == "1", f"score must be 1, got {parts[4]}"
            assert parts[5] == "+", f"orientation must be +, got {parts[5]}"
            # Strip "{sp}_" prefix — everything after the first underscore is busco_id
            busco_id = gene1_start.split("_", 1)[1]
            result[busco_id] = color
    return result


class TestGenerateLinksOneLinePerGene:
    """The .simple file must contain exactly one line per common BUSCO gene
    (no block aggregation) so that a gene's color can never be overridden by
    a block-level majority vote.
    """

    def test_one_line_per_common_gene(self, tmp_path):
        # 5 co-linear genes on the same chromosome pair — the old block
        # detector would have merged them into a single block.
        sp1 = {f"g{i}": _gene("chr1", i * 1000) for i in range(5)}
        sp2 = {f"g{i}": _gene("chrA", i * 1000) for i in range(5)}
        gene_colors = {f"g{i}": "#ff0000" for i in range(5)}
        path = str(tmp_path / "links.simple")

        generate_links_file(sp1, sp2, gene_colors, "sp1", "sp2", path)

        parsed = _parse_simple_file(path)
        assert set(parsed.keys()) == {f"g{i}" for i in range(5)}
        assert all(c == "#ff0000" for c in parsed.values())

    def test_mixed_colors_preserved_not_majority_voted(self, tmp_path):
        # A mixed run of 5 co-linear genes: 3 red, 2 blue. Old code would
        # have voted this a single red block. New code must preserve every
        # gene's individual color.
        sp1 = {f"g{i}": _gene("chr1", i * 1000) for i in range(5)}
        sp2 = {f"g{i}": _gene("chrA", i * 1000) for i in range(5)}
        gene_colors = {
            "g0": "#ff0000",
            "g1": "#ff0000",
            "g2": "#0000ff",
            "g3": "#ff0000",
            "g4": "#0000ff",
        }
        path = str(tmp_path / "links.simple")

        generate_links_file(sp1, sp2, gene_colors, "sp1", "sp2", path)

        parsed = _parse_simple_file(path)
        assert parsed == gene_colors

    def test_grey_written_before_colored(self, tmp_path):
        # JCVI renders later lines on top of earlier ones, so non-significant
        # (grey) lines must precede significant ones in the file.
        sp1 = {f"g{i}": _gene("chr1", i * 1000) for i in range(4)}
        sp2 = {f"g{i}": _gene("chrA", i * 1000) for i in range(4)}
        gene_colors = {
            "g0": "#ff0000",
            "g1": "lightgrey",
            "g2": "#00ff00",
            "g3": "lightgrey",
        }
        path = str(tmp_path / "links.simple")

        generate_links_file(sp1, sp2, gene_colors, "sp1", "sp2", path)

        with open(path) as f:
            lines = f.readlines()
        first_colored_idx = next(
            i for i, line in enumerate(lines) if not line.startswith("lightgrey*")
        )
        last_grey_idx = max(
            i for i, line in enumerate(lines) if line.startswith("lightgrey*")
        )
        assert last_grey_idx < first_colored_idx


class TestCrossFileGeneColorConsistency:
    """Regression: a gene shared between two species pairs must appear with
    the exact same color in both links.*.simple files.
    """

    def test_shared_gene_identical_color_across_pairs(self, tmp_path):
        # Shared ortholog `g_shared` present in all three species. The same
        # `gene_colors` mapping (derived from a single run-wide gene_to_chain)
        # must be used for both pairs, and the resulting lines must agree.
        sp1 = {"g_shared": _gene("chr1", 1000), "g_sp1_only": _gene("chr2", 2000)}
        sp2 = {"g_shared": _gene("chrA", 1000), "g_sp2_only": _gene("chrB", 2000)}
        sp3 = {"g_shared": _gene("chrX", 1000), "g_sp3_only": _gene("chrY", 2000)}
        gene_colors = {
            "g_shared": "#aabbcc",
            "g_sp1_only": "lightgrey",
            "g_sp2_only": "lightgrey",
            "g_sp3_only": "lightgrey",
        }

        path_ab = str(tmp_path / "links.sp1.sp2.simple")
        path_bc = str(tmp_path / "links.sp2.sp3.simple")
        generate_links_file(sp1, sp2, gene_colors, "sp1", "sp2", path_ab)
        generate_links_file(sp2, sp3, gene_colors, "sp2", "sp3", path_bc)

        parsed_ab = _parse_simple_file(path_ab)
        parsed_bc = _parse_simple_file(path_bc)
        assert parsed_ab["g_shared"] == parsed_bc["g_shared"] == "#aabbcc"

    def test_shared_gene_same_color_despite_different_neighbors(self, tmp_path):
        # g_shared is surrounded by red neighbors in the A-B pair and by
        # blue neighbors in the B-C pair. Under the old majority-vote code,
        # this was exactly the failure mode: a gene "inherited" its block's
        # dominant color, so it could appear red in A-B and blue in B-C.
        sp1 = {
            "g_shared": _gene("chr1", 2000),
            "g_red1": _gene("chr1", 1000),
            "g_red2": _gene("chr1", 3000),
        }
        sp2 = {
            "g_shared": _gene("chrA", 2000),
            "g_red1": _gene("chrA", 1000),
            "g_red2": _gene("chrA", 3000),
            "g_blue1": _gene("chrA", 4000),
            "g_blue2": _gene("chrA", 5000),
        }
        sp3 = {
            "g_shared": _gene("chrX", 3000),
            "g_blue1": _gene("chrX", 1000),
            "g_blue2": _gene("chrX", 2000),
        }
        gene_colors = {
            "g_shared": "#00ff00",  # green — a minority in both blocks
            "g_red1": "#ff0000",
            "g_red2": "#ff0000",
            "g_blue1": "#0000ff",
            "g_blue2": "#0000ff",
        }

        path_ab = str(tmp_path / "links.sp1.sp2.simple")
        path_bc = str(tmp_path / "links.sp2.sp3.simple")
        generate_links_file(sp1, sp2, gene_colors, "sp1", "sp2", path_ab)
        generate_links_file(sp2, sp3, gene_colors, "sp2", "sp3", path_bc)

        parsed_ab = _parse_simple_file(path_ab)
        parsed_bc = _parse_simple_file(path_bc)
        assert parsed_ab["g_shared"] == "#00ff00"
        assert parsed_bc["g_shared"] == "#00ff00"


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
            "save_chromosome_associations": {},
            "detect_algs_transitive": {"return_value": ([], [], {}, {}, [])},
            "apply_custom_colors": {"return_value": {"g1": "lightgrey"}},
            "apply_custom_colors_with_algs": {"return_value": {"g1": "lightgrey"}},
            "parse_custom_colors": {"return_value": {"g1": "#00ff00"}},
            "parse_custom_algs": {"return_value": {"g1": "ALG_A"}},
            "save_rearrangement_indices": {},
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

    def _call_run(self, tmp_path, custom_color_file="", skip_alg=False, **kwargs):
        run(
            assembly_busco="/fake/assembly_busco",
            assembly_fasta="/fake/assembly.fasta",
            busco_refs=["/fake/busco_reference_ref1"],
            accession_order="",
            manual_refs="/fake/ref1.fasta",
            output_dir=str(tmp_path / "output"),
            custom_color_file=custom_color_file,
            skip_alg=skip_alg,
            **kwargs,
        )

    def test_custom_with_skip_alg(self, tmp_path, run_mocks):
        self._call_run(tmp_path, custom_color_file="/fake/colors.tsv", skip_alg=True)

        run_mocks["apply_custom_colors"].assert_called()
        run_mocks["detect_algs_transitive"].assert_not_called()
        run_mocks["apply_custom_colors_with_algs"].assert_not_called()
        run_mocks["save_chromosome_associations"].assert_called_once()
        run_mocks["save_rearrangement_indices"].assert_called_once()

    def test_custom_without_skip_alg(self, tmp_path, run_mocks):
        self._call_run(tmp_path, custom_color_file="/fake/colors.tsv", skip_alg=False)

        run_mocks["apply_custom_colors_with_algs"].assert_called()
        run_mocks["detect_algs_transitive"].assert_called()
        run_mocks["apply_custom_colors"].assert_not_called()
        run_mocks["save_chromosome_associations"].assert_called_once()
        run_mocks["save_rearrangement_indices"].assert_called_once()

    def test_no_custom_colors(self, tmp_path, run_mocks):
        self._call_run(tmp_path, custom_color_file="", skip_alg=False)

        run_mocks["apply_custom_colors_with_algs"].assert_called()
        run_mocks["detect_algs_transitive"].assert_called()
        run_mocks["apply_custom_colors"].assert_not_called()
        run_mocks["save_chromosome_associations"].assert_called_once()
        run_mocks["save_rearrangement_indices"].assert_not_called()

    def test_alpha_passed_to_detect_algs_transitive(self, tmp_path, run_mocks):
        self._call_run(tmp_path, alpha=0.05)

        run_mocks["detect_algs_transitive"].assert_called_once()
        _, kwargs = run_mocks["detect_algs_transitive"].call_args
        assert kwargs["alpha"] == 0.05


def _assoc(sp1, sp2, chr1, chr2):
    return PairwiseAssociation(sp1, sp2, chr1, chr2, p_value=0.001, gene_count=10)


class TestEnumerateChains:
    def test_empty_input(self):
        assert enumerate_chains([], []) == []

    def test_two_species_single_pair(self):
        assocs = [_assoc("A", "B", "A1", "B1")]
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B1")}),
        ]
        chains = enumerate_chains(assocs, species_busco)
        assert chains == [[("A", "A1"), ("B", "B1")]]

    def test_two_species_multiple_pairs(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("A", "B", "A2", "B2"),
        ]
        species_busco = [
            ("A", {"g1": _gene("A1"), "g2": _gene("A2")}),
            ("B", {"g1": _gene("B1"), "g2": _gene("B2")}),
        ]
        chains = enumerate_chains(assocs, species_busco)
        assert len(chains) == 2
        assert [("A", "A1"), ("B", "B1")] in chains
        assert [("A", "A2"), ("B", "B2")] in chains

    def test_three_species_linear(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("B", "C", "B1", "C1"),
        ]
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B1")}),
            ("C", {"g1": _gene("C1")}),
        ]
        chains = enumerate_chains(assocs, species_busco)
        assert chains == [[("A", "A1"), ("B", "B1"), ("C", "C1")]]

    def test_three_species_downstream_branching(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("B", "C", "B1", "C2"),
            _assoc("B", "C", "B1", "C3"),
        ]
        species_busco = [
            ("A", {"g1": _gene("A1"), "g2": _gene("A1")}),
            ("B", {"g1": _gene("B1"), "g2": _gene("B1")}),
            ("C", {"g1": _gene("C2"), "g2": _gene("C3")}),
        ]
        chains = enumerate_chains(assocs, species_busco)
        assert len(chains) == 2
        assert [("A", "A1"), ("B", "B1"), ("C", "C2")] in chains
        assert [("A", "A1"), ("B", "B1"), ("C", "C3")] in chains

    def test_three_species_upstream_merging(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("A", "B", "A2", "B1"),
            _assoc("B", "C", "B1", "C2"),
        ]
        species_busco = [
            ("A", {"g1": _gene("A1"), "g2": _gene("A2")}),
            ("B", {"g1": _gene("B1"), "g2": _gene("B1")}),
            ("C", {"g1": _gene("C2"), "g2": _gene("C2")}),
        ]
        chains = enumerate_chains(assocs, species_busco)
        assert len(chains) == 2
        assert [("A", "A1"), ("B", "B1"), ("C", "C2")] in chains
        assert [("A", "A2"), ("B", "B1"), ("C", "C2")] in chains

    def test_three_species_dead_end(self):
        assocs = [_assoc("A", "B", "A1", "B1")]
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B1")}),
            ("C", {}),
        ]
        chains = enumerate_chains(assocs, species_busco)
        assert chains == [[("A", "A1"), ("B", "B1")]]

    def test_orphan_chain_later_species(self):
        assocs = [_assoc("B", "C", "B3", "C4")]
        species_busco = [
            ("A", {}),
            ("B", {"g1": _gene("B3")}),
            ("C", {"g1": _gene("C4")}),
        ]
        chains = enumerate_chains(assocs, species_busco)
        assert chains == [[("B", "B3"), ("C", "C4")]]

    def test_deterministic_ordering(self):
        assocs = [
            _assoc("A", "B", "A2", "B2"),
            _assoc("A", "B", "A1", "B1"),
        ]
        species_busco = [
            ("A", {"g1": _gene("A1"), "g2": _gene("A2")}),
            ("B", {"g1": _gene("B1"), "g2": _gene("B2")}),
        ]
        chains = enumerate_chains(assocs, species_busco)
        assert chains[0] == [("A", "A1"), ("B", "B1")]
        assert chains[1] == [("A", "A2"), ("B", "B2")]

    def test_mixed_dead_end_and_orphan(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("B", "C", "B3", "C4"),
        ]
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B1"), "g2": _gene("B3")}),
            ("C", {"g2": _gene("C4")}),
        ]
        chains = enumerate_chains(assocs, species_busco)
        assert len(chains) == 2
        assert [("A", "A1"), ("B", "B1")] in chains
        assert [("B", "B3"), ("C", "C4")] in chains

    def test_four_species_spanning_chain(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("B", "C", "B1", "C2"),
            _assoc("C", "D", "C2", "D3"),
        ]
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B1")}),
            ("C", {"g1": _gene("C2")}),
            ("D", {"g1": _gene("D3")}),
        ]
        chains = enumerate_chains(assocs, species_busco)
        assert chains == [[("A", "A1"), ("B", "B1"), ("C", "C2"), ("D", "D3")]]

    def test_diamond_converging_at_last_species(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("A", "B", "A1", "B2"),
            _assoc("B", "C", "B1", "C2"),
            _assoc("B", "C", "B2", "C2"),
        ]
        species_busco = [
            ("A", {"g1": _gene("A1"), "g2": _gene("A1")}),
            ("B", {"g1": _gene("B1"), "g2": _gene("B2")}),
            ("C", {"g1": _gene("C2"), "g2": _gene("C2")}),
        ]
        chains = enumerate_chains(assocs, species_busco)
        assert len(chains) == 2
        assert [("A", "A1"), ("B", "B1"), ("C", "C2")] in chains
        assert [("A", "A1"), ("B", "B2"), ("C", "C2")] in chains

    def test_subchain_removed_when_longer_chain_exists(self):
        """A gene breaking at a non-significant edge should not create a
        sub-chain that lets unrelated genes match."""
        assocs = [
            _assoc("A", "B", "A1", "B2"),
            _assoc("B", "C", "B2", "C4"),
        ]
        # g1 walks A1→B2 (sig) → C4 (sig) → full chain
        # g2 walks A1→B2 (sig) → C3 (NOT sig) → would emit sub-chain [(A,A1),(B,B2)]
        # That sub-chain is a prefix of the full chain and must be dropped.
        species_busco = [
            ("A", {"g1": _gene("A1"), "g2": _gene("A1")}),
            ("B", {"g1": _gene("B2"), "g2": _gene("B2")}),
            ("C", {"g1": _gene("C4"), "g2": _gene("C3")}),
        ]
        chains = enumerate_chains(assocs, species_busco)
        assert chains == [[("A", "A1"), ("B", "B2"), ("C", "C4")]]

    def test_cross_product_pruned_by_gene_evidence(self):
        """Shared node B2 should NOT produce cartesian product of chains."""
        assocs = [
            _assoc("A", "B", "A4", "B2"),
            _assoc("A", "B", "A9", "B2"),
            _assoc("B", "C", "B2", "C11"),
            _assoc("B", "C", "B2", "C14"),
        ]
        species_busco = [
            ("A", {"g1": _gene("A4"), "g2": _gene("A9")}),
            ("B", {"g1": _gene("B2"), "g2": _gene("B2")}),
            ("C", {"g1": _gene("C11"), "g2": _gene("C14")}),
        ]
        chains = enumerate_chains(assocs, species_busco)
        assert len(chains) == 2
        assert [("A", "A4"), ("B", "B2"), ("C", "C11")] in chains
        assert [("A", "A9"), ("B", "B2"), ("C", "C14")] in chains
        # Spurious chains must NOT exist:
        assert [("A", "A4"), ("B", "B2"), ("C", "C14")] not in chains
        assert [("A", "A9"), ("B", "B2"), ("C", "C11")] not in chains

    def test_chain_skips_absent_species(self):
        """Gene absent in middle species should not break the chain."""
        assocs = [_assoc("A", "C", "A1", "C1")]
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {}),
            ("C", {"g1": _gene("C1")}),
        ]
        chains = enumerate_chains(assocs, species_busco)
        assert chains == [[("A", "A1"), ("C", "C1")]]

    def test_skip_species_subchain_pruned_by_longer(self):
        """Chain skipping a species must be pruned when a longer chain
        covering that species exists."""
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("A", "C", "A1", "C1"),
            _assoc("B", "C", "B1", "C1"),
        ]
        # g1 present in all 3 → walks A1→B1→C1 (full chain)
        # g2 absent from B → walks A1→C1 (subsequence of full chain)
        species_busco = [
            ("A", {"g1": _gene("A1"), "g2": _gene("A1")}),
            ("B", {"g1": _gene("B1")}),
            ("C", {"g1": _gene("C1"), "g2": _gene("C1")}),
        ]
        chains = enumerate_chains(assocs, species_busco)
        # Within-chain order is not significant (consumers index by species),
        # so assert on the node set rather than the DFS canonicalization order.
        assert len(chains) == 1
        assert set(chains[0]) == {("A", "A1"), ("B", "B1"), ("C", "C1")}

    def test_order_independence(self):
        """Same species set, different input order → same chains."""
        # B has a fragmented genome for this ALG: gene g1 is NOT in B.
        # A-C and C-D have significant associations.
        assocs_order1 = [
            _assoc("A", "B", "A1", "B9"),
            _assoc("A", "C", "A1", "C1"),
            _assoc("A", "D", "A1", "D1"),
            _assoc("B", "C", "B9", "C1"),
            _assoc("B", "D", "B9", "D1"),
            _assoc("C", "D", "C1", "D1"),
        ]
        busco_order1 = [
            ("A", {"g1": _gene("A1")}),
            ("B", {}),
            ("C", {"g1": _gene("C1")}),
            ("D", {"g1": _gene("D1")}),
        ]

        # Same species, B moved to second position
        assocs_order2 = [
            _assoc("C", "B", "C1", "B9"),
            _assoc("C", "A", "C1", "A1"),
            _assoc("C", "D", "C1", "D1"),
            _assoc("B", "A", "B9", "A1"),
            _assoc("B", "D", "B9", "D1"),
            _assoc("A", "D", "A1", "D1"),
        ]
        busco_order2 = [
            ("C", {"g1": _gene("C1")}),
            ("B", {}),
            ("A", {"g1": _gene("A1")}),
            ("D", {"g1": _gene("D1")}),
        ]

        chains1 = enumerate_chains(assocs_order1, busco_order1)
        chains2 = enumerate_chains(assocs_order2, busco_order2)

        # Normalize: sort each chain's nodes by species name, then sort chains
        def normalize(chains):
            return sorted(sorted(c, key=lambda n: n[0]) for c in chains)

        assert normalize(chains1) == normalize(chains2)
        assert len(chains1) == 1

    def test_order_independence_partial_edges(self):
        """Order independence with partial pairwise edges (not all pairs significant)."""
        # A-B and B-C significant, but A-C NOT significant.
        # Gene g1 present in A, B, C → should form chain [A,B,C].
        assocs_order1 = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("B", "C", "B1", "C1"),
        ]
        busco_order1 = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B1")}),
            ("C", {"g1": _gene("C1")}),
        ]

        # Reversed species order
        assocs_order2 = [
            _assoc("C", "B", "C1", "B1"),
            _assoc("B", "A", "B1", "A1"),
        ]
        busco_order2 = [
            ("C", {"g1": _gene("C1")}),
            ("B", {"g1": _gene("B1")}),
            ("A", {"g1": _gene("A1")}),
        ]

        chains1 = enumerate_chains(assocs_order1, busco_order1)
        chains2 = enumerate_chains(assocs_order2, busco_order2)

        def normalize(chains):
            return sorted(sorted(c, key=lambda n: n[0]) for c in chains)

        assert normalize(chains1) == normalize(chains2)
        assert len(chains1) == 1

    def test_order_independence_four_species_sparse(self):
        """Order independence with 4 species, only consecutive edges significant."""
        assocs_fwd = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("B", "C", "B1", "C1"),
            _assoc("C", "D", "C1", "D1"),
        ]
        busco_fwd = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B1")}),
            ("C", {"g1": _gene("C1")}),
            ("D", {"g1": _gene("D1")}),
        ]

        # Scrambled order: C, A, D, B
        # Same 3 significant pairs, species1/species2 follow i<j in new order
        assocs_scrambled = [
            _assoc("C", "B", "C1", "B1"),
            _assoc("C", "D", "C1", "D1"),
            _assoc("A", "B", "A1", "B1"),
        ]
        busco_scrambled = [
            ("C", {"g1": _gene("C1")}),
            ("A", {"g1": _gene("A1")}),
            ("D", {"g1": _gene("D1")}),
            ("B", {"g1": _gene("B1")}),
        ]

        chains_fwd = enumerate_chains(assocs_fwd, busco_fwd)
        chains_scr = enumerate_chains(assocs_scrambled, busco_scrambled)

        def normalize(chains):
            return sorted(sorted(c, key=lambda n: n[0]) for c in chains)

        assert normalize(chains_fwd) == normalize(chains_scr)
        norm = normalize(chains_fwd)
        assert len(norm) == 1
        assert norm[0] == [("A", "A1"), ("B", "B1"), ("C", "C1"), ("D", "D1")]

    def test_order_independence_detect_algs_transitive(self):
        """Full pipeline: detect_algs_transitive gives same chains regardless
        of species_busco ordering."""
        busco_order1 = [
            ("A", {"g1": _gene("A1"), "g2": _gene("A2")}),
            ("B", {"g1": _gene("B1"), "g2": _gene("B2")}),
            ("C", {"g1": _gene("C1"), "g2": _gene("C2")}),
        ]
        busco_order2 = [
            ("C", {"g1": _gene("C1"), "g2": _gene("C2")}),
            ("A", {"g1": _gene("A1"), "g2": _gene("A2")}),
            ("B", {"g1": _gene("B1"), "g2": _gene("B2")}),
        ]

        _, chains1, mapping1, _, _ = detect_algs_transitive(
            busco_order1, min_genes=1, min_chain_genes=0
        )
        _, chains2, mapping2, _, _ = detect_algs_transitive(
            busco_order2, min_genes=1, min_chain_genes=0
        )

        def normalize(chains):
            return sorted(sorted(c, key=lambda n: n[0]) for c in chains)

        assert normalize(chains1) == normalize(chains2)
        # Gene-to-chain mapping should match (same gene → same chain content)
        for gene_id in mapping1:
            if mapping1[gene_id] == -1:
                assert mapping2[gene_id] == -1
            else:
                chain1_nodes = set(
                    (sp, ch) for sp, ch in chains1[mapping1[gene_id]]
                )
                chain2_nodes = set(
                    (sp, ch) for sp, ch in chains2[mapping2[gene_id]]
                )
                assert chain1_nodes == chain2_nodes, f"gene {gene_id} mapped differently"


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
    @patch(f"{STATS_MODULE}.detect_algs_pairwise_raw")
    def test_two_species_returns_gene_to_chain(self, mock_pairwise):
        sp1 = {"g1": _gene("A1"), "g2": _gene("A2")}
        sp2 = {"g1": _gene("B1"), "g2": _gene("B2")}
        mock_pairwise.return_value = (
            [
                _assoc("sp1", "sp2", "A1", "B1"),
                _assoc("sp1", "sp2", "A2", "B2"),
            ],
            [],
        )

        assocs, chains, gene_to_chain, chain_colors, _ = detect_algs_transitive(
            [("sp1", sp1), ("sp2", sp2)], min_chain_genes=1
        )

        assert len(assocs) == 2
        assert len(chains) == 2
        assert isinstance(gene_to_chain, dict)
        assert all(isinstance(v, int) for v in gene_to_chain.values())
        assert gene_to_chain["g1"] != gene_to_chain["g2"]
        assert gene_to_chain["g1"] in chain_colors
        assert gene_to_chain["g2"] in chain_colors

    @patch(f"{STATS_MODULE}.detect_algs_pairwise_raw")
    def test_each_chain_gets_palette_color(self, mock_pairwise):
        sp1 = {"g1": _gene("A1"), "g2": _gene("A2")}
        sp2 = {"g1": _gene("B1"), "g2": _gene("B2")}
        mock_pairwise.return_value = (
            [
                _assoc("sp1", "sp2", "A1", "B1"),
                _assoc("sp1", "sp2", "A2", "B2"),
            ],
            [],
        )

        _, _, _, chain_colors, _ = detect_algs_transitive([("sp1", sp1), ("sp2", sp2)], min_chain_genes=1)

        for chain_id, color in chain_colors.items():
            assert color == ALG_PALETTE[chain_id % len(ALG_PALETTE)]

    @patch(f"{STATS_MODULE}.detect_algs_pairwise_raw")
    def test_three_species_branching_distinct_chains(self, mock_pairwise):
        sp1 = {"g1": _gene("A1"), "g2": _gene("A1")}
        sp2 = {"g1": _gene("B1"), "g2": _gene("B1")}
        sp3 = {"g1": _gene("C1"), "g2": _gene("C2")}
        mock_pairwise.side_effect = [
            ([_assoc("sp1", "sp2", "A1", "B1")], []),
            ([], []),  # sp1-sp3 (non-adjacent, no effect on chains)
            (
                [
                    _assoc("sp2", "sp3", "B1", "C1"),
                    _assoc("sp2", "sp3", "B1", "C2"),
                ],
                [],
            ),
        ]

        _, _, gene_to_chain, chain_colors, _ = detect_algs_transitive(
            [("sp1", sp1), ("sp2", sp2), ("sp3", sp3)], min_chain_genes=1
        )

        assert gene_to_chain["g1"] != gene_to_chain["g2"]
        assert len(chain_colors) == 2

    @patch(f"{STATS_MODULE}.detect_algs_pairwise_raw")
    def test_no_associations_returns_empty(self, mock_pairwise):
        sp1 = {"g1": _gene("A1")}
        sp2 = {"g1": _gene("B1")}
        mock_pairwise.return_value = ([], [])

        assocs, chains, gene_to_chain, chain_colors, chr_assocs = (
            detect_algs_transitive([("sp1", sp1), ("sp2", sp2)], min_chain_genes=1)
        )

        assert assocs == []
        assert chains == []
        assert chain_colors == {}
        assert all(v == -1 for v in gene_to_chain.values())
        assert chr_assocs == []

    @patch(f"{STATS_MODULE}.detect_algs_pairwise_raw")
    def test_single_species_returns_empty(self, mock_pairwise):
        sp1 = {"g1": _gene("A1")}

        assocs, chains, gene_to_chain, chain_colors, chr_assocs = (
            detect_algs_transitive([("sp1", sp1)])
        )

        mock_pairwise.assert_not_called()
        assert assocs == []
        assert chains == []
        assert gene_to_chain == {}
        assert chain_colors == {}
        assert chr_assocs == []


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
            (
                [_assoc("sp1", "sp2", "A1", "B1"), _assoc("sp1", "sp2", "A1", "B2")],
                [],
            ),
            ([], []),  # sp1-sp3 (non-adjacent, no effect on chains)
            (
                [_assoc("sp2", "sp3", "B1", "C2"), _assoc("sp2", "sp3", "B2", "C3")],
                [],
            ),
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
            ([_assoc("sp1", "sp2", "A1", "B1")], []),
            ([], []),  # sp1-sp3 (non-adjacent, no effect on chains)
            ([_assoc("sp2", "sp3", "B1", "C2")], []),
        ]

    @patch(f"{STATS_MODULE}.detect_algs_pairwise_raw")
    def test_branching_genes_get_different_colors(
        self, mock_pairwise, branching_species, branching_pairwise
    ):
        mock_pairwise.side_effect = branching_pairwise

        _, chains, gene_to_chain, chain_colors, _ = detect_algs_transitive(
            branching_species, min_chain_genes=1
        )

        assert gene_to_chain["g_b1"] != gene_to_chain["g_b2"]
        assert gene_to_chain["g_b1"] >= 0
        assert gene_to_chain["g_b2"] >= 0
        assert gene_to_chain["g_nosig"] == -1

        sp1_busco = branching_species[0][1]
        sp2_busco = branching_species[1][1]
        colors_ab = apply_custom_colors_with_algs(
            sp1_busco, sp2_busco, gene_to_chain, chain_colors, {},
            chains=chains, sp1_name="sp1", sp2_name="sp2",
        )
        assert colors_ab["g_b1"] != colors_ab["g_b2"]
        assert colors_ab["g_b1"] != "lightgrey"
        assert colors_ab["g_b2"] != "lightgrey"
        assert colors_ab["g_nosig"] == "lightgrey"

    @patch(f"{STATS_MODULE}.detect_algs_pairwise_raw")
    def test_linear_chain_same_color_across_pairs(
        self, mock_pairwise, linear_species, linear_pairwise
    ):
        mock_pairwise.side_effect = linear_pairwise

        _, chains, gene_to_chain, chain_colors, _ = detect_algs_transitive(
            linear_species, min_chain_genes=1
        )

        colors_ab = apply_custom_colors_with_algs(
            linear_species[0][1], linear_species[1][1], gene_to_chain, chain_colors, {},
            chains=chains, sp1_name="sp1", sp2_name="sp2",
        )
        colors_bc = apply_custom_colors_with_algs(
            linear_species[1][1], linear_species[2][1], gene_to_chain, chain_colors, {},
            chains=chains, sp1_name="sp2", sp2_name="sp3",
        )

        assert colors_ab["g_sig"] == colors_bc["g_sig"]
        assert colors_ab["g_sig"] != "lightgrey"

    @patch(f"{STATS_MODULE}.detect_algs_pairwise_raw")
    def test_linear_chain_with_custom_colors(
        self, mock_pairwise, linear_species, linear_pairwise
    ):
        mock_pairwise.side_effect = linear_pairwise

        _, chains, gene_to_chain, chain_colors, _ = detect_algs_transitive(
            linear_species, min_chain_genes=1
        )

        custom = {"g_sig": "#00ff00"}
        sp1_busco = linear_species[0][1]
        sp2_busco = linear_species[1][1]

        colors = apply_custom_colors_with_algs(
            sp1_busco, sp2_busco, gene_to_chain, chain_colors, custom,
            chains=chains, sp1_name="sp1", sp2_name="sp2",
        )

        assert colors["g_sig"] == "#00ff00"
        assert colors["g_nosig"] == "lightgrey"


class TestTwoSpeciesBackwardCompat:
    """Verify 2-species chain-based pipeline produces correct coloring."""

    @patch(f"{STATS_MODULE}.detect_algs_pairwise_raw")
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
        mock_pairwise.return_value = (
            [
                _assoc("sp1", "sp2", "A1", "B1"),
                _assoc("sp1", "sp2", "A2", "B2"),
            ],
            [],
        )

        _, chains, gene_to_chain, chain_colors, _ = detect_algs_transitive(
            [("sp1", sp1), ("sp2", sp2)], min_chain_genes=1
        )

        colors = apply_custom_colors_with_algs(
            sp1, sp2, gene_to_chain, chain_colors, {},
            chains=chains, sp1_name="sp1", sp2_name="sp2",
        )

        assert colors["g1"] == colors["g3"]
        assert colors["g1"] != colors["g2"]
        assert colors["g1"] in ALG_PALETTE
        assert colors["g2"] in ALG_PALETTE
        assert colors["g_nosig"] == "lightgrey"
        assert len(chains) == 2
        for cid, color in chain_colors.items():
            assert color == ALG_PALETTE[cid % len(ALG_PALETTE)]
