"""Tests for coloring functions and run-wide color exports."""

from csv import DictReader
from unittest.mock import patch

import pytest

from hobrac.jcvi_synteny.chains import build_gene_chain_mapping, enumerate_chains
from hobrac.jcvi_synteny.coloring import (
    apply_custom_colors,
    apply_custom_colors_with_algs,
)
from hobrac.jcvi_synteny.models import (
    ALG_PALETTE,
    AlgDetection,
    BuscoGene,
    PairwiseAssociation,
)
from hobrac.jcvi_synteny.output import generate_links_file, run
from hobrac.jcvi_synteny.statistics import detect_algs_transitive

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
        display = {"g1": 0, "g2": -1}
        chain_colors = {0: "#ff0000"}
        chains = [[("sp1", "chr1"), ("sp2", "chrA")]]
        return sp1, sp2, display, chain_colors, chains

    def test_empty_custom_colors_gives_chain_colors(self, chain_fixture):
        sp1, sp2, display, chain_colors, chains = chain_fixture
        result = apply_custom_colors_with_algs(
            sp1,
            sp2,
            display,
            chain_colors,
            {},
            chains=chains,
            sp1_name="sp1",
            sp2_name="sp2",
        )
        assert result["g1"] == chain_colors[0]

    def test_all_genes_in_custom(self, chain_fixture):
        sp1, sp2, display, chain_colors, chains = chain_fixture
        result = apply_custom_colors_with_algs(
            sp1,
            sp2,
            display,
            chain_colors,
            {"g1": "#00ff00", "g2": "#0000ff"},
            chains=chains,
            sp1_name="sp1",
            sp2_name="sp2",
        )
        # Custom colors are gated by chain membership: g1 is on a chain and
        # listed → custom; g2 is off-chain → lightgrey even though listed.
        assert result["g1"] == "#00ff00"
        assert result["g2"] == "lightgrey"

    def test_no_genes_in_custom(self, chain_fixture):
        sp1, sp2, display, chain_colors, chains = chain_fixture
        result = apply_custom_colors_with_algs(
            sp1,
            sp2,
            display,
            chain_colors,
            {},
            chains=chains,
            sp1_name="sp1",
            sp2_name="sp2",
        )
        assert result["g1"] == "#ff0000"
        assert result["g2"] == "lightgrey"

    def test_only_selected_display_chain_can_color_a_species_pair(self):
        sp_b = {"g1": _gene("B1")}
        sp_c = {"g1": _gene("C1")}
        sp_d = {"g1": _gene("D1")}
        chains = [
            [("A", "A1"), ("B", "B1"), ("C", "C1")],
            [("C", "C1"), ("D", "D1")],
        ]
        display = {"g1": 0}
        chain_colors = {0: "#ff0000", 1: "#0000ff"}
        custom = {"g1": "#00ff00"}

        colors_bc = apply_custom_colors_with_algs(
            sp_b,
            sp_c,
            display,
            chain_colors,
            custom,
            chains=chains,
            sp1_name="B",
            sp2_name="C",
        )
        colors_cd = apply_custom_colors_with_algs(
            sp_c,
            sp_d,
            display,
            chain_colors,
            custom,
            chains=chains,
            sp1_name="C",
            sp2_name="D",
        )

        assert colors_bc["g1"] == "#00ff00"
        assert colors_cd["g1"] == "lightgrey"

    def test_chain_covers_pair_canonical_order(self):
        """Chain in canonical (lex) order must still cover a pair whose
        species_busco order is non-lex."""
        sp_z = {"g1": _gene("Z1")}
        sp_a = {"g1": _gene("A1")}
        # Chain in canonical order: A before Z
        chains = [[("A", "A1"), ("Z", "Z1")]]
        display = {"g1": 0}
        chain_colors = {0: "#ff0000"}

        # species_busco order is Z, A (non-lex)
        colors = apply_custom_colors_with_algs(
            sp_z,
            sp_a,
            display,
            chain_colors,
            {},
            chains=chains,
            sp1_name="Z",
            sp2_name="A",
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
            assert parts[1] == gene1_start, (
                "start1 != end1 (blocks must be single-gene)"
            )
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
        # Shared ortholog `g_shared` is present in all three species. The same
        # run-wide selected display chain color is used for both pairs, so the
        # resulting lines must agree.
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


class TestRunColorExports:
    @pytest.fixture
    def run_inputs(self, tmp_path):
        # Two significant chains, plus one gene crossing between them.
        shared_rows = "".join(
            f"g{i}\tComplete\tchr{i // 10 + 1}\t{i * 100}\t{i * 100 + 90}\n"
            for i in range(20)
        )
        assembly_busco = tmp_path / "assembly.tsv"
        assembly_busco.write_text(
            shared_rows + "unlinked\tComplete\tchr1\t2000\t2090\n"
        )
        ref_dir = tmp_path / "busco_reference_ref1"
        ref_run = ref_dir / "run_busco"
        ref_run.mkdir(parents=True)
        (ref_run / "full_table.tsv").write_text(
            shared_rows + "unlinked\tComplete\tchr2\t2000\t2090\n"
        )
        assembly_fasta = tmp_path / "assembly.fasta"
        assembly_fasta.write_text(f">chr1\n{'A' * 2200}\n>chr2\n{'A' * 2200}\n")
        colors = tmp_path / "colors.tsv"
        colors.write_text(
            "g0\t#123456\tALG_A\ng1\t#abcdef\tALG_B\nunlinked\t#fedcba\tALG_C\n"
        )
        return {
            "assembly_busco": str(assembly_busco),
            "assembly_fasta": str(assembly_fasta),
            "busco_refs": [str(ref_dir)],
            "accession_order": "",
            "manual_refs": "ref1.fasta",
            "output_dir": str(tmp_path / "output"),
        }, str(colors)

    def _read_tsv(self, path):
        with open(path) as f:
            return list(DictReader(f, delimiter="\t"))

    def test_custom_colors_do_not_label_inferred_chains(self, run_inputs):
        inputs, colors = run_inputs
        output = run(**inputs, custom_color_file=colors)

        algs = self._read_tsv(output["algs"])
        assert {
            (row["assembly"], row["ref1"], row["n_genes"], row["color"]) for row in algs
        } == {("chr1", "chr1", "10", "-"), ("chr2", "chr2", "10", "-")}

        genes = {row["gene"]: row for row in self._read_tsv(output["gene_chains"])}
        assert genes["g0"]["chain_ids"] == genes["g1"]["chain_ids"] != "-"
        assert genes["g0"]["display_chain_id"] == genes["g1"]["display_chain_id"] != "-"
        assert genes["g0"]["custom_alg_id"] == "ALG_A"
        assert genes["g1"]["custom_alg_id"] == "ALG_B"
        assert genes["g0"]["color"] == "#123456"
        assert genes["g1"]["color"] == "#abcdef"
        assert genes["g2"]["color"] == "lightgrey"
        assert genes["unlinked"]["color"] == "lightgrey"
        assert _parse_simple_file(output["links_files"][0]) == {
            gene: row["color"] for gene, row in genes.items()
        }

    def test_no_custom_colors_exports_plot_palette(self, run_inputs):
        inputs, _ = run_inputs
        output = run(**inputs)

        algs = self._read_tsv(output["algs"])
        palette = {row["chain_id"]: row["color"] for row in algs}
        assert set(palette.values()) == set(ALG_PALETTE[:2])
        genes = self._read_tsv(output["gene_chains"])
        assert list(genes[0])[:5] == [
            "chain_ids",
            "custom_alg_id",
            "gene",
            "color",
            "display_chain_id",
        ]
        for row in genes:
            assert row["color"] == palette.get(row["display_chain_id"], "lightgrey")
        assert _parse_simple_file(output["links_files"][0]) == {
            row["gene"]: row["color"] for row in genes
        }

    def test_skip_alg_keeps_custom_colors_without_inferred_chains(self, run_inputs):
        inputs, colors = run_inputs
        output = run(**inputs, custom_color_file=colors, skip_alg=True)

        assert self._read_tsv(output["algs"]) == []
        genes = {row["gene"]: row for row in self._read_tsv(output["gene_chains"])}
        assert {row["chain_ids"] for row in genes.values()} == {"-"}
        assert {row["display_chain_id"] for row in genes.values()} == {"-"}
        assert genes["g0"]["custom_alg_id"] == "ALG_A"
        assert genes["g1"]["custom_alg_id"] == "ALG_B"
        assert genes["unlinked"]["custom_alg_id"] == "ALG_C"
        assert genes["g0"]["color"] == "#123456"
        assert genes["g1"]["color"] == "#abcdef"
        assert genes["g2"]["color"] == "lightgrey"
        assert genes["unlinked"]["color"] == "#fedcba"
        assert _parse_simple_file(output["links_files"][0]) == {
            gene: row["color"] for gene, row in genes.items()
        }


def _assoc(sp1, sp2, chr1, chr2):
    return PairwiseAssociation(sp1, sp2, chr1, chr2, p_value=0.001, gene_count=10)


def _chain(*nodes):
    return tuple(sorted(nodes))


class TestEnumerateChains:
    def test_empty_input(self):
        assert enumerate_chains([], []) == {}

    def test_observed_genes_are_the_chain_support(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("B", "C", "B1", "C1"),
        ]
        species_busco = [
            ("A", {"g1": _gene("A1"), "g2": _gene("A1")}),
            ("B", {"g1": _gene("B1"), "g2": _gene("B1")}),
            ("C", {"g1": _gene("C1"), "g2": _gene("C1")}),
        ]

        assert enumerate_chains(assocs, species_busco) == {
            _chain(("A", "A1"), ("B", "B1"), ("C", "C1")): {"g1", "g2"}
        }

    def test_star_retains_every_longest_node_set_with_distinct_supporters(self):
        assocs = [
            _assoc("B", "A", "B1", "A1"),
            _assoc("B", "C", "B1", "C1"),
            _assoc("B", "D", "B1", "D1"),
        ]
        species_busco = [
            ("A", {"all": _gene("A1"), "missing_d": _gene("A1")}),
            ("B", {"all": _gene("B1"), "missing_d": _gene("B1")}),
            ("C", {"all": _gene("C1"), "missing_d": _gene("C1")}),
            ("D", {"all": _gene("D1")}),
        ]

        assert enumerate_chains(assocs, species_busco) == {
            _chain(("A", "A1"), ("B", "B1"), ("C", "C1")): {
                "all",
                "missing_d",
            },
            _chain(("A", "A1"), ("B", "B1"), ("D", "D1")): {"all"},
            _chain(("B", "B1"), ("C", "C1"), ("D", "D1")): {"all"},
        }

    def test_gene_votes_once_when_multiple_traversals_have_same_node_set(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("B", "C", "B1", "C1"),
            _assoc("C", "D", "C1", "D1"),
            _assoc("D", "A", "D1", "A1"),
        ]
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B1")}),
            ("C", {"g1": _gene("C1")}),
            ("D", {"g1": _gene("D1")}),
        ]

        assert enumerate_chains(assocs, species_busco) == {
            _chain(
                ("A", "A1"),
                ("B", "B1"),
                ("C", "C1"),
                ("D", "D1"),
            ): {"g1"}
        }

    @patch("hobrac.jcvi_synteny.chains._MAX_DFS_VISITS", 1)
    def test_exact_search_bound_fails_instead_of_returning_partial_evidence(self):
        with pytest.raises(RuntimeError):
            enumerate_chains(
                [_assoc("A", "B", "A1", "B1")],
                [
                    ("A", {"g1": _gene("A1")}),
                    ("B", {"g1": _gene("B1")}),
                ],
            )

    def test_candidate_identity_and_support_are_input_order_and_name_invariant(self):
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("B", "C", "B1", "C1"),
            _assoc("A", "B", "A2", "B2"),
        ]
        original = [
            ("A", {"g_full": _gene("A1"), "g_pair": _gene("A2")}),
            ("B", {"g_full": _gene("B1"), "g_pair": _gene("B2")}),
            ("C", {"g_full": _gene("C1")}),
        ]
        renamed_and_reordered = [
            ("C", {"renamed_full": _gene("C1")}),
            (
                "B",
                {
                    "renamed_pair": _gene("B2"),
                    "renamed_full": _gene("B1"),
                },
            ),
            (
                "A",
                {
                    "renamed_pair": _gene("A2"),
                    "renamed_full": _gene("A1"),
                },
            ),
        ]
        reversed_assocs = [
            _assoc(a.species2, a.species1, a.chr2, a.chr1) for a in reversed(assocs)
        ]

        first = enumerate_chains(assocs, original)
        second = enumerate_chains(reversed_assocs, renamed_and_reordered)

        assert list(first) == list(second)
        assert [len(supporters) for supporters in first.values()] == [
            len(supporters) for supporters in second.values()
        ]

    def test_support_threshold_at_five_runs_before_subset_pruning(self):
        short = [f"short{i}" for i in range(5)]
        long = [f"long{i}" for i in range(4)]
        species_busco = [
            (
                "A",
                {gene: _gene("A1") for gene in short + long},
            ),
            (
                "B",
                {gene: _gene("B1") for gene in short + long},
            ),
            (
                "C",
                {gene: _gene("C1") for gene in long},
            ),
        ]
        assocs = [
            _assoc("A", "B", "A1", "B1"),
            _assoc("B", "C", "B1", "C1"),
        ]

        assert enumerate_chains(assocs, species_busco) == {
            _chain(("A", "A1"), ("B", "B1"), ("C", "C1")): set(long)
        }
        assert enumerate_chains(
            assocs,
            species_busco,
            min_chain_genes=5,
        ) == {_chain(("A", "A1"), ("B", "B1")): set(short)}


class TestBuildGeneChainMapping:
    def test_returns_all_candidates_and_selected_display_chain(self):
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B1")}),
            ("C", {"g1": _gene("C1")}),
        ]
        chains = [
            [("A", "A1"), ("B", "B1")],
            [("A", "A1"), ("B", "B1"), ("C", "C1")],
        ]

        candidates, display = build_gene_chain_mapping(
            species_busco, chains, {0: 3, 1: 7}
        )

        assert candidates == {"g1": [0, 1]}
        assert display == {"g1": 1}

    def test_missing_species_is_compatible_but_adds_no_visible_link(self):
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B1")}),
            ("X", {}),
            ("C", {"g1": _gene("C1")}),
            ("D", {"g1": _gene("D1")}),
        ]
        chains = [
            [("A", "A1"), ("B", "B1"), ("X", "X1"), ("C", "C1")],
            [("A", "A1"), ("B", "B1"), ("C", "C1"), ("D", "D1")],
        ]

        candidates, display = build_gene_chain_mapping(
            species_busco, chains, {0: 4, 1: 4}
        )

        assert candidates["g1"] == [0, 1]
        assert display["g1"] == 1

    def test_observed_support_precedes_visible_coverage(self):
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B1")}),
            ("C", {"g1": _gene("C1")}),
        ]
        chains = [
            [("A", "A1"), ("B", "B1")],
            [("A", "A1"), ("B", "B1"), ("C", "C1")],
        ]

        candidates, display = build_gene_chain_mapping(
            species_busco, chains, {0: 8, 1: 7}
        )

        assert candidates["g1"] == [0, 1]
        assert display["g1"] == 0

    def test_chain_id_is_the_last_display_tiebreaker(self):
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B1")}),
        ]
        chains = [
            [("A", "A1"), ("B", "B1")],
            [("A", "A1"), ("B", "B1")],
        ]

        candidates, display = build_gene_chain_mapping(
            species_busco, chains, {0: 5, 1: 5}
        )

        assert candidates["g1"] == [0, 1]
        assert display["g1"] == 0

    def test_assignments_never_reinforce_later_display_choices(self):
        species_busco = [
            (
                "A",
                {
                    "ambiguous_a": _gene("A1"),
                    "ambiguous_b": _gene("A1"),
                    "exact": _gene("A1"),
                },
            ),
            (
                "B",
                {
                    "ambiguous_a": _gene("B1"),
                    "ambiguous_b": _gene("B1"),
                    "exact": _gene("B1"),
                },
            ),
            ("C", {"exact": _gene("C3")}),
        ]
        chains = [
            [("A", "A1"), ("B", "B1"), ("C", "C2")],
            [("A", "A1"), ("B", "B1"), ("C", "C3")],
        ]

        candidates, display = build_gene_chain_mapping(
            species_busco, chains, {0: 5, 1: 5}
        )

        assert candidates == {
            "ambiguous_a": [0, 1],
            "ambiguous_b": [0, 1],
            "exact": [1],
        }
        assert display == {
            "ambiguous_a": 0,
            "ambiguous_b": 0,
            "exact": 1,
        }

    def test_unmatched_and_empty_chains_keep_eligible_genes_explicit(self):
        species_busco = [
            ("A", {"g1": _gene("A1")}),
            ("B", {"g1": _gene("B9")}),
        ]

        candidates, display = build_gene_chain_mapping(
            species_busco,
            [[("A", "A1"), ("B", "B1")]],
            {0: 1},
        )
        empty_candidates, empty_display = build_gene_chain_mapping(
            species_busco, [], {}
        )

        assert candidates == {"g1": []}
        assert display == {"g1": -1}
        assert empty_candidates == {"g1": []}
        assert empty_display == {"g1": -1}

    def test_gene_in_single_species_is_excluded(self):
        species_busco = [
            ("A", {"shared": _gene("A1"), "lone": _gene("A9")}),
            ("B", {"shared": _gene("B1")}),
        ]

        candidates, display = build_gene_chain_mapping(
            species_busco,
            [[("A", "A1"), ("B", "B1")]],
            {0: 1},
        )

        assert "lone" not in candidates
        assert "lone" not in display


class TestDetectAlgsTransitive:
    @patch(f"{STATS_MODULE}.detect_algs_pairwise_raw")
    def test_observed_support_is_frozen_before_ambiguous_display_selection(
        self, mock_pairwise
    ):
        sp1 = {
            "g1": _gene("A1"),
            "g2": _gene("A1"),
            "g3": _gene("A1"),
            "ambiguous": _gene("A1"),
        }
        sp2 = {
            "g1": _gene("B1"),
            "g2": _gene("B1"),
            "g3": _gene("B1"),
            "ambiguous": _gene("B1"),
        }
        sp3 = {
            "g1": _gene("C2"),
            "g2": _gene("C2"),
            "g3": _gene("C3"),
        }
        mock_pairwise.side_effect = [
            ([_assoc("sp1", "sp2", "A1", "B1")], []),
            ([], []),
            (
                [
                    _assoc("sp2", "sp3", "B1", "C2"),
                    _assoc("sp2", "sp3", "B1", "C3"),
                ],
                [],
            ),
        ]

        result = detect_algs_transitive(
            [("sp1", sp1), ("sp2", sp2), ("sp3", sp3)],
            min_chain_genes=1,
        )

        assert isinstance(result, AlgDetection)
        assert result.chain_support == {0: 2, 1: 1}
        assert result.gene_to_chains == {
            "ambiguous": [0, 1],
            "g1": [0],
            "g2": [0],
            "g3": [1],
        }
        assert result.gene_to_display_chain["ambiguous"] == 0
        assert result.chain_colors == {
            0: ALG_PALETTE[0],
            1: ALG_PALETTE[1],
        }
        assert len(result.associations) == 3

    @patch(f"{STATS_MODULE}.detect_algs_pairwise_raw")
    def test_no_associations_returns_empty_dataclass_fields(self, mock_pairwise):
        sp1 = {"g1": _gene("A1")}
        sp2 = {"g1": _gene("B1")}
        mock_pairwise.return_value = ([], [])

        result = detect_algs_transitive([("sp1", sp1), ("sp2", sp2)], min_chain_genes=1)

        assert isinstance(result, AlgDetection)
        assert result.associations == []
        assert result.chains == []
        assert result.chain_support == {}
        assert result.gene_to_chains == {"g1": []}
        assert result.gene_to_display_chain == {"g1": -1}
        assert result.chain_colors == {}
        assert result.chromosome_associations == []

    @patch(f"{STATS_MODULE}.detect_algs_pairwise_raw")
    def test_single_species_does_not_run_pairwise_detection(self, mock_pairwise):
        result = detect_algs_transitive([("sp1", {"g1": _gene("A1")})])

        mock_pairwise.assert_not_called()
        assert isinstance(result, AlgDetection)
        assert result.associations == []
        assert result.chains == []
        assert result.chain_support == {}
        assert result.gene_to_chains == {}
        assert result.gene_to_display_chain == {}
        assert result.chain_colors == {}
        assert result.chromosome_associations == []


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

        detection = detect_algs_transitive(branching_species, min_chain_genes=1)

        display = detection.gene_to_display_chain
        assert display["g_b1"] != display["g_b2"]
        assert display["g_b1"] >= 0
        assert display["g_b2"] >= 0
        assert detection.gene_to_chains["g_nosig"] == []
        assert display["g_nosig"] == -1

        sp1_busco = branching_species[0][1]
        sp2_busco = branching_species[1][1]
        colors_ab = apply_custom_colors_with_algs(
            sp1_busco,
            sp2_busco,
            display,
            detection.chain_colors,
            {},
            chains=detection.chains,
            sp1_name="sp1",
            sp2_name="sp2",
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

        detection = detect_algs_transitive(linear_species, min_chain_genes=1)

        colors_ab = apply_custom_colors_with_algs(
            linear_species[0][1],
            linear_species[1][1],
            detection.gene_to_display_chain,
            detection.chain_colors,
            {},
            chains=detection.chains,
            sp1_name="sp1",
            sp2_name="sp2",
        )
        colors_bc = apply_custom_colors_with_algs(
            linear_species[1][1],
            linear_species[2][1],
            detection.gene_to_display_chain,
            detection.chain_colors,
            {},
            chains=detection.chains,
            sp1_name="sp2",
            sp2_name="sp3",
        )

        assert colors_ab["g_sig"] == colors_bc["g_sig"]
        assert colors_ab["g_sig"] != "lightgrey"

    @patch(f"{STATS_MODULE}.detect_algs_pairwise_raw")
    def test_linear_chain_with_custom_colors(
        self, mock_pairwise, linear_species, linear_pairwise
    ):
        mock_pairwise.side_effect = linear_pairwise

        detection = detect_algs_transitive(linear_species, min_chain_genes=1)

        custom = {"g_sig": "#00ff00"}
        sp1_busco = linear_species[0][1]
        sp2_busco = linear_species[1][1]

        colors = apply_custom_colors_with_algs(
            sp1_busco,
            sp2_busco,
            detection.gene_to_display_chain,
            detection.chain_colors,
            custom,
            chains=detection.chains,
            sp1_name="sp1",
            sp2_name="sp2",
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

        detection = detect_algs_transitive(
            [("sp1", sp1), ("sp2", sp2)], min_chain_genes=1
        )

        colors = apply_custom_colors_with_algs(
            sp1,
            sp2,
            detection.gene_to_display_chain,
            detection.chain_colors,
            {},
            chains=detection.chains,
            sp1_name="sp1",
            sp2_name="sp2",
        )

        assert colors["g1"] == colors["g3"]
        assert colors["g1"] != colors["g2"]
        assert colors["g1"] in ALG_PALETTE
        assert colors["g2"] in ALG_PALETTE
        assert colors["g_nosig"] == "lightgrey"
        assert len(detection.chains) == 2
        for cid, color in detection.chain_colors.items():
            assert color == ALG_PALETTE[cid % len(ALG_PALETTE)]
