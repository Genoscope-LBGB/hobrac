"""Deduplication of NCBI Datasets rows into the precompute_mash input list."""

from hobrac.dedup_ncbi import deduplicate


def test_distinct_assemblies_sharing_a_name_are_all_kept():
    rows = [
        ("Version 1", "GCA_049487705.1", "", "41880"),
        ("Version 1", "GCA_049488355.1", "", "41886"),
        ("Version 1", "GCF_049487705.1", "GCA_049487705.1", "41880"),
    ]
    assert deduplicate(rows) == [
        ("GCA_049487705.1", "Version 1", "41880"),
        ("GCA_049488355.1", "Version 1", "41886"),
    ]


def test_gcf_row_collapses_into_paired_gca_even_when_names_differ():
    rows = [
        ("mLepTim1.pri", "GCF_033115175.1", "GCA_033115175.1", "62621"),
        ("mLepEur2.pri", "GCA_033115175.1", "GCF_033115175.1", "62621"),
    ]
    assert deduplicate(rows) == [("GCA_033115175.1", "mLepEur2.pri", "62621")]


def test_old_gca_version_paired_with_gcf_yields_to_newest_gca():
    rows = [
        ("asm.1", "GCA_943734635.1", "GCF_943734635.1", "7"),
        ("asm.1", "GCF_943734635.1", "GCA_943734635.1", "7"),
        ("asm.2", "GCA_943734635.2", "", "7"),
        ("other", "GCA_000000001.1", "", "8"),
    ]
    assert deduplicate(rows) == [
        ("GCA_943734635.2", "asm.2", "7"),
        ("GCA_000000001.1", "other", "8"),
    ]


def test_gcf_without_genbank_counterpart_is_kept():
    rows = [("ASM14184v1", "GCF_000141845.1", "", "412030")]
    assert deduplicate(rows) == [("GCF_000141845.1", "ASM14184v1", "412030")]
