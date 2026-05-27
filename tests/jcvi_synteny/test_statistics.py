"""Tests for detect_algs_pairwise_raw: corrected p-values and significance."""

from hobrac.jcvi_synteny.models import BuscoGene
from hobrac.jcvi_synteny.statistics import detect_algs_pairwise_raw


def _gene(chromosome, start=0):
    return BuscoGene(busco_id="", chromosome=chromosome, start=start, end=start + 100)


def _build_busco(gene_assignments):
    return {
        gid: _gene(chrom, i * 100) for i, (gid, chrom) in enumerate(gene_assignments)
    }


class TestDetectAlgsPairwiseRaw:
    def test_returns_both_significant_and_non_significant(self):
        # chr3→chrA: 50 genes all on the same pair → very significant.
        sp1_genes = [(f"s{i}", "chr3") for i in range(50)]
        sp2_genes = [(f"s{i}", "chrA") for i in range(50)]

        # chr1 and chr2: 10 genes each, split evenly between chrA and chrB
        # → no enrichment, not significant after Bonferroni correction.
        sp1_genes += [(f"a{i}", "chr1") for i in range(10)]
        sp2_genes += [(f"a{i}", "chrA") for i in range(5)]
        sp2_genes += [(f"a{i}", "chrB") for i in range(5, 10)]

        sp1_genes += [(f"b{i}", "chr2") for i in range(10)]
        sp2_genes += [(f"b{i}", "chrA") for i in range(5)]
        sp2_genes += [(f"b{i}", "chrB") for i in range(5, 10)]

        sp1 = _build_busco(sp1_genes)
        sp2 = _build_busco(sp2_genes)

        significant, all_tested = detect_algs_pairwise_raw(
            sp1, sp2, "sp1", "sp2", alpha=0.01, min_genes=5
        )

        assert len(all_tested) > 0

        sig_flags = [a.significant for a in all_tested]
        assert True in sig_flags, "expected at least one significant pair"
        assert False in sig_flags, "expected at least one non-significant pair"

        sig_chrs = {(a.chr1, a.chr2) for a in significant}
        for a in all_tested:
            if a.significant:
                assert (a.chr1, a.chr2) in sig_chrs
            else:
                assert (a.chr1, a.chr2) not in sig_chrs

    def test_corrected_p_value_is_bonferroni(self):
        # Create data with multiple testable pairs so num_tests > 1
        sp1_genes = [(f"g{i}", "chr1") for i in range(10)]
        sp2_genes = [(f"g{i}", "chrA") for i in range(10)]
        sp1_genes += [(f"h{i}", "chr2") for i in range(8)]
        sp2_genes += [(f"h{i}", "chrB") for i in range(8)]

        sp1 = _build_busco(sp1_genes)
        sp2 = _build_busco(sp2_genes)

        _, all_tested = detect_algs_pairwise_raw(sp1, sp2, "sp1", "sp2", min_genes=5)

        num_tests = len(all_tested)
        assert num_tests >= 2, "need multiple testable pairs for this test"

        for a in all_tested:
            expected = min(a.p_value * num_tests, 1.0)
            assert a.corrected_p_value == expected, (
                f"corrected_p_value {a.corrected_p_value} != "
                f"min({a.p_value} * {num_tests}, 1.0) = {expected}"
            )

    def test_corrected_p_value_capped_at_one(self):
        # Force a pair with a high raw p-value so correction would exceed 1.0
        # chr1-chrA: 5 genes (testable but weak association)
        # chr1-chrB: 5 genes (testable, creates noise)
        sp1_genes = [(f"a{i}", "chr1") for i in range(10)]
        sp2_genes = [(f"a{i}", "chrA") for i in range(5)]
        sp2_genes += [(f"a{i}", "chrB") for i in range(5, 10)]

        sp1 = _build_busco(sp1_genes)
        sp2 = _build_busco(sp2_genes)

        _, all_tested = detect_algs_pairwise_raw(sp1, sp2, "sp1", "sp2", min_genes=5)

        for a in all_tested:
            assert a.corrected_p_value <= 1.0

    def test_significance_matches_bonferroni_threshold(self):
        # Strong association: all genes on same chromosome pair
        sp1_genes = [(f"g{i}", "chr1") for i in range(15)]
        sp2_genes = [(f"g{i}", "chrA") for i in range(15)]
        # Weak pair to create a second test
        sp1_genes += [(f"h{i}", "chr2") for i in range(5)]
        sp2_genes += [(f"h{i}", "chrB") for i in range(5)]

        sp1 = _build_busco(sp1_genes)
        sp2 = _build_busco(sp2_genes)

        _, all_tested = detect_algs_pairwise_raw(
            sp1, sp2, "sp1", "sp2", alpha=0.01, min_genes=5
        )

        num_tests = len(all_tested)
        threshold = 0.01 / num_tests

        for a in all_tested:
            expected_sig = a.p_value <= threshold
            assert a.significant == expected_sig, (
                f"({a.chr1}, {a.chr2}): significant={a.significant} but "
                f"p_value={a.p_value} vs threshold={threshold}"
            )

    def test_species_names_propagated(self):
        sp1 = _build_busco([(f"g{i}", "chr1") for i in range(10)])
        sp2 = _build_busco([(f"g{i}", "chrA") for i in range(10)])

        _, all_tested = detect_algs_pairwise_raw(
            sp1, sp2, "Homo_sapiens", "Mus_musculus", min_genes=5
        )

        assert len(all_tested) > 0
        for a in all_tested:
            assert a.species1 == "Homo_sapiens"
            assert a.species2 == "Mus_musculus"

    def test_empty_input_returns_empty(self):
        significant, all_tested = detect_algs_pairwise_raw({}, {}, "sp1", "sp2")
        assert significant == []
        assert all_tested == []

    def test_no_common_genes_returns_empty(self):
        sp1 = _build_busco([("g1", "chr1")])
        sp2 = _build_busco([("g2", "chrA")])

        significant, all_tested = detect_algs_pairwise_raw(sp1, sp2, "sp1", "sp2")
        assert significant == []
        assert all_tested == []

    def test_below_min_genes_returns_empty(self):
        # 4 genes on same pair, but min_genes=5 → not testable
        sp1 = _build_busco([(f"g{i}", "chr1") for i in range(4)])
        sp2 = _build_busco([(f"g{i}", "chrA") for i in range(4)])

        significant, all_tested = detect_algs_pairwise_raw(
            sp1, sp2, "sp1", "sp2", min_genes=5
        )
        assert significant == []
        assert all_tested == []
