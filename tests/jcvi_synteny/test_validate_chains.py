"""Tests for validate_chains: pruning chain nodes by significant link count."""

from hobrac.jcvi_synteny.chains import validate_chains
from hobrac.jcvi_synteny.models import BuscoGene, ChromosomeAssociation


def _sig(sp1, chr1, sp2, chr2):
    return ChromosomeAssociation(
        species1=sp1,
        species2=sp2,
        chr1=chr1,
        chr2=chr2,
        p_value=0.001,
        corrected_p_value=0.001,
        gene_count=10,
        significant=True,
    )


def _nonsig(sp1, chr1, sp2, chr2):
    return ChromosomeAssociation(
        species1=sp1,
        species2=sp2,
        chr1=chr1,
        chr2=chr2,
        p_value=0.5,
        corrected_p_value=1.0,
        gene_count=3,
        significant=False,
    )


class TestValidateChains:
    def test_chain_length_2_always_passes(self):
        chain = [("sp1", "chr1"), ("sp2", "chrA")]
        assocs = [_sig("sp1", "chr1", "sp2", "chrA")]
        result = validate_chains([chain], assocs)
        assert result == [chain]

    def test_chain_length_3_passes_with_consecutive_only(self):
        chain = [("sp1", "chr1"), ("sp2", "chrA"), ("sp3", "chrX")]
        assocs = [
            _sig("sp1", "chr1", "sp2", "chrA"),
            _sig("sp2", "chrA", "sp3", "chrX"),
        ]
        result = validate_chains([chain], assocs)
        assert result == [chain]

    def test_chain_length_4_strict_prunes_weak_node(self):
        chain = [("sp1", "c1"), ("sp2", "c2"), ("sp3", "c3"), ("sp4", "c4")]
        assocs = [
            _sig("sp1", "c1", "sp2", "c2"),
            _sig("sp2", "c2", "sp3", "c3"),
            _sig("sp3", "c3", "sp4", "c4"),
            _sig("sp1", "c1", "sp3", "c3"),
            _sig("sp2", "c2", "sp4", "c4"),
            # sp1-c1 has links to sp2-c2 and sp3-c3 → 2 links, passes (threshold=2)
            # sp2-c2 has links to sp1-c1, sp3-c3, sp4-c4 → 3 links, passes
            # sp3-c3 has links to sp2-c2, sp4-c4, sp1-c1 → 3 links, passes
            # sp4-c4 has links to sp3-c3, sp2-c2 → 2 links, passes
        ]
        result = validate_chains([chain], assocs)
        assert result == [chain]

    def test_chain_length_4_strict_drops_node_with_one_link(self):
        chain = [("sp1", "c1"), ("sp2", "c2"), ("sp3", "c3"), ("sp4", "c4")]
        assocs = [
            _sig("sp1", "c1", "sp2", "c2"),
            _sig("sp2", "c2", "sp3", "c3"),
            _sig("sp3", "c3", "sp4", "c4"),
            _sig("sp1", "c1", "sp3", "c3"),
            # Round 1 (n=4, threshold=2):
            # sp4-c4: 1 link (sp3-c3) → pruned
            # sp1-c1: 2 links → passes; sp2-c2: 2 links → passes
            # sp3-c3: 3 links (sp2-c2, sp1-c1, sp4-c4) → passes
            # Round 2 (n=3, threshold=1): all pass
        ]
        result = validate_chains([chain], assocs)
        assert result == [[("sp1", "c1"), ("sp2", "c2"), ("sp3", "c3")]]

    def test_chain_length_4_permissive_keeps_node_with_one_link(self):
        chain = [("sp1", "c1"), ("sp2", "c2"), ("sp3", "c3"), ("sp4", "c4")]
        assocs = [
            _sig("sp1", "c1", "sp2", "c2"),
            _sig("sp2", "c2", "sp3", "c3"),
            _sig("sp3", "c3", "sp4", "c4"),
        ]
        result = validate_chains([chain], assocs, permissive=True)
        assert result == [chain]

    def test_prune_splits_chain_into_two(self):
        chain = [
            ("sp1", "c1"),
            ("sp2", "c2"),
            ("sp3", "c3"),
            ("sp4", "c4"),
            ("sp5", "c5"),
        ]
        # Middle node sp3-c3 has no significant links at all → pruned, splits chain
        assocs = [
            _sig("sp1", "c1", "sp2", "c2"),
            _sig("sp4", "c4", "sp5", "c5"),
            _sig("sp1", "c1", "sp4", "c4"),
            _sig("sp1", "c1", "sp5", "c5"),
            _sig("sp2", "c2", "sp4", "c4"),
            _sig("sp2", "c2", "sp5", "c5"),
            _sig("sp4", "c4", "sp5", "c5"),
        ]
        result = validate_chains([chain], assocs, permissive=True)
        assert [("sp1", "c1"), ("sp2", "c2")] in result
        assert [("sp4", "c4"), ("sp5", "c5")] in result
        assert len(result) == 2

    def test_no_significant_links_drops_chain(self):
        chain = [("sp1", "c1"), ("sp2", "c2"), ("sp3", "c3")]
        assocs = [_nonsig("sp1", "c1", "sp2", "c2")]
        result = validate_chains([chain], assocs)
        assert result == []

    def test_empty_chains(self):
        assert validate_chains([], []) == []

    def test_empty_associations_drops_all(self):
        chain = [("sp1", "c1"), ("sp2", "c2")]
        result = validate_chains([chain], [])
        assert result == []

    def test_nonsignificant_associations_not_counted(self):
        chain = [("sp1", "c1"), ("sp2", "c2"), ("sp3", "c3"), ("sp4", "c4")]
        assocs = [
            _sig("sp1", "c1", "sp2", "c2"),
            _nonsig("sp1", "c1", "sp3", "c3"),
            _nonsig("sp1", "c1", "sp4", "c4"),
            _sig("sp2", "c2", "sp3", "c3"),
            _sig("sp3", "c3", "sp4", "c4"),
            _sig("sp2", "c2", "sp4", "c4"),
        ]
        # sp1-c1: only 1 sig link (sp2-c2), threshold=2 → pruned
        # remaining: sp2, sp3, sp4 each have 2 sig links → kept
        result = validate_chains([chain], assocs)
        assert result == [[("sp2", "c2"), ("sp3", "c3"), ("sp4", "c4")]]

    def test_single_node_after_prune_dropped(self):
        chain = [("sp1", "c1"), ("sp2", "c2"), ("sp3", "c3")]
        # sp3-c3 has no links → pruned; sp1-c1 and sp2-c2 share a link → kept
        assocs = [_sig("sp2", "c2", "sp1", "c1")]
        result = validate_chains([chain], assocs, permissive=True)
        assert result == [[("sp1", "c1"), ("sp2", "c2")]]

    def test_multiple_chains_validated_independently(self):
        chain1 = [("sp1", "c1"), ("sp2", "c2")]
        chain2 = [("sp1", "c3"), ("sp2", "c4")]
        assocs = [
            _sig("sp1", "c1", "sp2", "c2"),
            # no sig link for chain2
        ]
        result = validate_chains([chain1, chain2], assocs)
        assert result == [chain1]

    def test_iterative_pruning_removes_node_surviving_on_pruned_links(self):
        # Chain of 6: A and B each link only to D, E, F (not to each other).
        # D, E, F each link only to A and B. First pass (threshold=3):
        # A count=3 (D,E,F), B count=3 (D,E,F), D/E/F count=2 → pruned.
        # Second pass on [A,B] (threshold=1): A has 0 links to B → pruned.
        chain = [
            ("sp1", "c1"),
            ("sp2", "c2"),
            ("sp3", "c3"),
            ("sp4", "c4"),
            ("sp5", "c5"),
            ("sp6", "c6"),
        ]
        assocs = [
            _sig("sp1", "c1", "sp4", "c4"),
            _sig("sp1", "c1", "sp5", "c5"),
            _sig("sp1", "c1", "sp6", "c6"),
            _sig("sp2", "c2", "sp4", "c4"),
            _sig("sp2", "c2", "sp5", "c5"),
            _sig("sp2", "c2", "sp6", "c6"),
        ]
        result = validate_chains([chain], assocs)
        assert result == []

    def test_duplicate_subchains_deduplicated(self):
        chain1 = [("sp1", "c1"), ("sp2", "c2"), ("sp3", "c3")]
        chain2 = [("sp1", "c1"), ("sp2", "c2"), ("sp3", "c4")]
        assocs = [
            _sig("sp1", "c1", "sp2", "c2"),
            # sp3-c3 and sp3-c4 have no cross-species links → pruned from both
        ]
        # Both chains pruned to [("sp1","c1"),("sp2","c2")] — should appear once
        result = validate_chains([chain1, chain2], assocs, permissive=True)
        assert result == [[("sp1", "c1"), ("sp2", "c2")]]

    def test_permissive_does_not_bypass_min_chain_genes(self):
        """Permissive mode relaxes link threshold but min_chain_genes still applies."""
        chain = [
            ("sp1", "c1"),
            ("sp2", "c2"),
            ("sp3", "c3"),
            ("sp4", "c4"),
            ("sp5", "c5"),
        ]
        # All consecutive pairs significant — passes permissive validation (threshold=1)
        assocs = [
            _sig("sp1", "c1", "sp2", "c2"),
            _sig("sp2", "c2", "sp3", "c3"),
            _sig("sp3", "c3", "sp4", "c4"),
            _sig("sp4", "c4", "sp5", "c5"),
        ]

        def _busco(chrom):
            return BuscoGene(chromosome=chrom, start=0, end=100, busco_id="x")

        # Only 3 genes support this chain — below default min_chain_genes=5
        species_busco = [
            ("sp1", {"g1": _busco("c1"), "g2": _busco("c1"), "g3": _busco("c1")}),
            ("sp2", {"g1": _busco("c2"), "g2": _busco("c2"), "g3": _busco("c2")}),
            ("sp3", {"g1": _busco("c3"), "g2": _busco("c3"), "g3": _busco("c3")}),
            ("sp4", {"g1": _busco("c4"), "g2": _busco("c4"), "g3": _busco("c4")}),
            ("sp5", {"g1": _busco("c5"), "g2": _busco("c5"), "g3": _busco("c5")}),
        ]

        # Without min_chain_genes, permissive keeps chain intact
        result_no_min = validate_chains(
            [chain], assocs, permissive=True, min_chain_genes=0,
        )
        assert result_no_min == [chain]

        # With min_chain_genes=5, chain dropped despite permissive=True
        result_with_min = validate_chains(
            [chain],
            assocs,
            permissive=True,
            min_chain_genes=5,
            species_busco=species_busco,
        )
        assert result_with_min == []

    def test_subchain_pruned_regardless_of_internal_order(self):
        """A shorter chain whose nodes are a subset of a longer chain must be
        pruned even when their internal orderings differ (canonical path
        direction may flip)."""
        chain_long = [("sp1", "c1"), ("sp2", "c2"), ("sp3", "c3"), ("sp4", "c4")]
        chain_short = [("sp3", "c3"), ("sp2", "c2"), ("sp1", "c1")]
        assocs = [
            _sig("sp1", "c1", "sp2", "c2"),
            _sig("sp2", "c2", "sp3", "c3"),
            _sig("sp3", "c3", "sp4", "c4"),
            _sig("sp1", "c1", "sp3", "c3"),
            _sig("sp2", "c2", "sp4", "c4"),
            _sig("sp1", "c1", "sp4", "c4"),
        ]
        result = validate_chains([chain_long, chain_short], assocs)
        assert result == [chain_long]

    def test_duplicate_chains_different_order_deduplicated(self):
        """Two chains with identical node sets but different orderings must be
        deduplicated to a single chain."""
        chain_a = [("sp1", "c1"), ("sp2", "c2"), ("sp3", "c3")]
        chain_b = [("sp3", "c3"), ("sp2", "c2"), ("sp1", "c1")]
        assocs = [
            _sig("sp1", "c1", "sp2", "c2"),
            _sig("sp2", "c2", "sp3", "c3"),
            _sig("sp1", "c1", "sp3", "c3"),
        ]
        result = validate_chains([chain_a, chain_b], assocs)
        assert len(result) == 1
        assert set(result[0]) == {("sp1", "c1"), ("sp2", "c2"), ("sp3", "c3")}

    def test_validation_pruning_creates_subchain_then_removed(self):
        """When validation prunes a node from a long chain, creating a shorter
        chain that is a subset of another surviving chain, the shorter one
        must be removed."""
        chain_full = [("sp1", "c1"), ("sp2", "c2"), ("sp3", "c3"), ("sp4", "c4")]
        chain_partial = [("sp1", "c1"), ("sp2", "c2"), ("sp3", "c3"), ("sp4", "c5")]
        assocs = [
            _sig("sp1", "c1", "sp2", "c2"),
            _sig("sp2", "c2", "sp3", "c3"),
            _sig("sp3", "c3", "sp4", "c4"),
            _sig("sp1", "c1", "sp3", "c3"),
            _sig("sp2", "c2", "sp4", "c4"),
            _sig("sp1", "c1", "sp4", "c4"),
            # sp4-c5 has no links → pruned from chain_partial
            # chain_partial becomes [sp1-c1, sp2-c2, sp3-c3] ⊂ chain_full
        ]
        result = validate_chains([chain_full, chain_partial], assocs)
        assert result == [chain_full]

    def test_min_chain_genes_filters_weak_subchains(self):
        # Two chains: one with strong gene support, one with weak
        chain_strong = [("sp1", "c1"), ("sp2", "c2")]
        chain_weak = [("sp1", "c3"), ("sp2", "c4")]
        assocs = [
            _sig("sp1", "c1", "sp2", "c2"),
            _sig("sp1", "c3", "sp2", "c4"),
        ]

        def _busco(chrom):
            return BuscoGene(chromosome=chrom, start=0, end=100, busco_id="x")

        species_busco = [
            ("sp1", {"g1": _busco("c1"), "g2": _busco("c1"), "g3": _busco("c3")}),
            ("sp2", {"g1": _busco("c2"), "g2": _busco("c2"), "g3": _busco("c4")}),
        ]
        # g1, g2 match chain_strong (2 genes); g3 matches chain_weak (1 gene)
        result = validate_chains(
            [chain_strong, chain_weak],
            assocs,
            permissive=True,
            min_chain_genes=2,
            species_busco=species_busco,
        )
        assert result == [chain_strong]
