"""Tests for coloring functions and run() branching logic."""

from unittest.mock import patch

import pytest

from hobrac.jcvi_synteny import (
    BuscoGene,
    PairwiseAssociation,
    apply_custom_colors,
    apply_custom_colors_with_algs,
    build_gene_colors_from_algs,
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
    def alg_fixture(self):
        sp1 = {"g1": _gene("chr1"), "g2": _gene("chr3")}
        sp2 = {"g1": _gene("chrA"), "g2": _gene("chrC")}
        chr_to_alg = {("sp1", "chr1"): 0, ("sp2", "chrA"): 0}
        alg_colors = {0: "#ff0000"}
        sig = [PairwiseAssociation("sp1", "sp2", "chr1", "chrA", 0.001, 10)]
        return sp1, sp2, chr_to_alg, alg_colors, sig

    def test_empty_custom_colors_matches_build_gene_colors(self, alg_fixture):
        sp1, sp2, chr_to_alg, alg_colors, sig = alg_fixture
        hybrid = apply_custom_colors_with_algs(
            sp1, sp2, "sp1", "sp2", chr_to_alg, alg_colors, sig, {}
        )
        alg_only = build_gene_colors_from_algs(
            sp1, sp2, "sp1", "sp2", chr_to_alg, alg_colors, sig
        )
        assert hybrid == alg_only

    def test_all_genes_in_custom(self, alg_fixture):
        sp1, sp2, chr_to_alg, alg_colors, sig = alg_fixture
        result = apply_custom_colors_with_algs(
            sp1,
            sp2,
            "sp1",
            "sp2",
            chr_to_alg,
            alg_colors,
            sig,
            {"g1": "#00ff00", "g2": "#0000ff"},
        )
        assert result["g1"] == "#00ff00"
        assert result["g2"] == "lightgrey"

    def test_no_genes_in_custom(self, alg_fixture):
        sp1, sp2, chr_to_alg, alg_colors, sig = alg_fixture
        result = apply_custom_colors_with_algs(
            sp1, sp2, "sp1", "sp2", chr_to_alg, alg_colors, sig, {}
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
            "detect_algs_transitive": {"return_value": ([], {}, {})},
            "build_alg_association_list": {"return_value": []},
            "apply_custom_colors": {"return_value": {"g1": "lightgrey"}},
            "apply_custom_colors_with_algs": {"return_value": {"g1": "lightgrey"}},
            "build_gene_colors_from_algs": {"return_value": {"g1": "lightgrey"}},
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
        run_mocks["build_gene_colors_from_algs"].assert_not_called()

    def test_custom_without_skip_alg(self, tmp_path, run_mocks):
        self._call_run(tmp_path, custom_color_file="/fake/colors.tsv", skip_alg=False)

        run_mocks["apply_custom_colors_with_algs"].assert_called()
        run_mocks["detect_algs_transitive"].assert_called()
        run_mocks["apply_custom_colors"].assert_not_called()
        run_mocks["build_gene_colors_from_algs"].assert_not_called()

    def test_no_custom_colors(self, tmp_path, run_mocks):
        self._call_run(tmp_path, custom_color_file="", skip_alg=False)

        run_mocks["build_gene_colors_from_algs"].assert_called()
        run_mocks["detect_algs_transitive"].assert_called()
        run_mocks["apply_custom_colors"].assert_not_called()
        run_mocks["apply_custom_colors_with_algs"].assert_not_called()
