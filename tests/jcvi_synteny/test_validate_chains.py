"""Behavioral tests for graph validation of observed chain evidence."""

from hobrac.jcvi_synteny.chains import validate_chains
from hobrac.jcvi_synteny.models import ChromosomeAssociation


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


def _chain(*nodes):
    return tuple(sorted(nodes))


class TestValidateChains:
    def test_valid_chain_preserves_observed_supporters(self):
        chain = _chain(("sp1", "c1"), ("sp2", "c2"), ("sp3", "c3"))
        evidence = {chain: {"g1", "g2"}}
        assocs = [
            _sig("sp1", "c1", "sp2", "c2"),
            _sig("sp2", "c2", "sp3", "c3"),
        ]

        assert validate_chains(evidence, assocs) == evidence

    def test_strict_degree_rule_prunes_a_weak_node(self):
        original = _chain(
            ("sp1", "c1"),
            ("sp2", "c2"),
            ("sp3", "c3"),
            ("sp4", "c4"),
        )
        evidence = {original: {"observed"}}
        assocs = [
            _sig("sp1", "c1", "sp2", "c2"),
            _sig("sp2", "c2", "sp3", "c3"),
            _sig("sp3", "c3", "sp4", "c4"),
            _sig("sp1", "c1", "sp3", "c3"),
        ]

        assert validate_chains(evidence, assocs) == {
            _chain(("sp1", "c1"), ("sp2", "c2"), ("sp3", "c3")): {"observed"}
        }

    def test_permissive_mode_keeps_nodes_with_one_significant_link(self):
        chain = _chain(
            ("sp1", "c1"),
            ("sp2", "c2"),
            ("sp3", "c3"),
            ("sp4", "c4"),
        )
        evidence = {chain: {"g1"}}
        assocs = [
            _sig("sp1", "c1", "sp2", "c2"),
            _sig("sp2", "c2", "sp3", "c3"),
            _sig("sp3", "c3", "sp4", "c4"),
        ]

        assert validate_chains(evidence, assocs, permissive=True) == evidence

    def test_connected_survivors_stay_together_across_input_list_gaps(self):
        original = _chain(
            ("sp1", "c1"),
            ("sp2", "unused"),
            ("sp3", "c3"),
        )
        evidence = {original: {"g1"}}

        assert validate_chains(
            evidence,
            [_sig("sp1", "c1", "sp3", "c3")],
            permissive=True,
        ) == {_chain(("sp1", "c1"), ("sp3", "c3")): {"g1"}}

    def test_disconnected_significant_pairs_split_by_graph_connectivity(self):
        original = _chain(
            ("sp1", "c1"),
            ("sp2", "c2"),
            ("sp3", "c3"),
            ("sp4", "c4"),
        )
        evidence = {original: {"g1"}}

        assert validate_chains(
            evidence,
            [
                _sig("sp1", "c1", "sp2", "c2"),
                _sig("sp3", "c3", "sp4", "c4"),
            ],
            permissive=True,
        ) == {
            _chain(("sp1", "c1"), ("sp2", "c2")): {"g1"},
            _chain(("sp3", "c3"), ("sp4", "c4")): {"g1"},
        }

    def test_all_longest_paths_from_a_surviving_star_are_retained(self):
        original = _chain(
            ("A", "A1"),
            ("B", "B1"),
            ("C", "C1"),
            ("D", "D1"),
        )
        evidence = {original: {"g1"}}
        assocs = [
            _sig("B", "B1", "A", "A1"),
            _sig("B", "B1", "C", "C1"),
            _sig("B", "B1", "D", "D1"),
        ]

        assert validate_chains(evidence, assocs, permissive=True) == {
            _chain(("A", "A1"), ("B", "B1"), ("C", "C1")): {"g1"},
            _chain(("A", "A1"), ("B", "B1"), ("D", "D1")): {"g1"},
            _chain(("B", "B1"), ("C", "C1"), ("D", "D1")): {"g1"},
        }

    def test_split_candidates_union_distinct_original_supporters(self):
        shared = (("sp1", "c1"), ("sp2", "c2"))
        evidence = {
            _chain(*shared, ("sp3", "c3")): {"g1", "g2", "shared"},
            _chain(*shared, ("sp4", "c4")): {"g3", "g4", "shared"},
        }

        assert validate_chains(
            evidence,
            [_sig("sp1", "c1", "sp2", "c2")],
            permissive=True,
            min_chain_genes=5,
        ) == {_chain(*shared): {"g1", "g2", "g3", "g4", "shared"}}

    def test_missing_or_partially_matching_genes_never_appear_as_support(self):
        chain = _chain(("sp1", "c1"), ("sp2", "c2"))
        evidence = {chain: {"observed"}}

        assert (
            validate_chains(
                evidence,
                [_sig("sp1", "c1", "sp2", "c2")],
            )
            == evidence
        )
        assert (
            validate_chains(
                evidence,
                [_sig("sp1", "c1", "sp2", "c2")],
                min_chain_genes=2,
            )
            == {}
        )

    def test_support_threshold_is_inclusive_at_five(self):
        keep = _chain(("sp1", "c1"), ("sp2", "c2"))
        drop = _chain(("sp1", "c3"), ("sp2", "c4"))
        evidence = {
            keep: {f"g{i}" for i in range(5)},
            drop: {f"h{i}" for i in range(4)},
        }

        assert validate_chains(
            evidence,
            [
                _sig("sp1", "c1", "sp2", "c2"),
                _sig("sp1", "c3", "sp2", "c4"),
            ],
            min_chain_genes=5,
        ) == {keep: {f"g{i}" for i in range(5)}}

    def test_nonsignificant_links_do_not_contribute_to_degree(self):
        original = _chain(
            ("sp1", "c1"),
            ("sp2", "c2"),
            ("sp3", "c3"),
            ("sp4", "c4"),
        )
        evidence = {original: {"g1"}}
        assocs = [
            _sig("sp1", "c1", "sp2", "c2"),
            _nonsig("sp1", "c1", "sp3", "c3"),
            _nonsig("sp1", "c1", "sp4", "c4"),
            _sig("sp2", "c2", "sp3", "c3"),
            _sig("sp3", "c3", "sp4", "c4"),
            _sig("sp2", "c2", "sp4", "c4"),
        ]

        assert validate_chains(evidence, assocs) == {
            _chain(("sp2", "c2"), ("sp3", "c3"), ("sp4", "c4")): {"g1"}
        }

    def test_pruning_is_iterative(self):
        original = _chain(
            ("sp1", "c1"),
            ("sp2", "c2"),
            ("sp3", "c3"),
            ("sp4", "c4"),
            ("sp5", "c5"),
            ("sp6", "c6"),
        )
        evidence = {original: {"g1"}}
        assocs = [
            _sig("sp1", "c1", "sp4", "c4"),
            _sig("sp1", "c1", "sp5", "c5"),
            _sig("sp1", "c1", "sp6", "c6"),
            _sig("sp2", "c2", "sp4", "c4"),
            _sig("sp2", "c2", "sp5", "c5"),
            _sig("sp2", "c2", "sp6", "c6"),
        ]

        assert validate_chains(evidence, assocs) == {}

    def test_empty_evidence_or_no_significant_links_returns_empty(self):
        chain = _chain(("sp1", "c1"), ("sp2", "c2"))

        assert validate_chains({}, []) == {}
        assert validate_chains({chain: {"g1"}}, []) == {}
        assert (
            validate_chains(
                {chain: {"g1"}},
                [_nonsig("sp1", "c1", "sp2", "c2")],
            )
            == {}
        )
