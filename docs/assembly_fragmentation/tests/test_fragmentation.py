"""Offline tests for level construction, sampling and report metrics."""

import gzip
import tempfile
import unittest
from pathlib import Path

from ..scripts import fragmentation


def seq_row(
    accession, role, chr_name="", unit="Primary Assembly", location="Chromosome"
):
    return {
        "genbank_accession": accession,
        "role": role,
        "chr_name": chr_name,
        "assembly_unit": unit,
        "assigned_molecule_location_type": location,
    }


def read_levels(directory):
    result = {}
    for level in fragmentation.LEVELS:
        with gzip.open(directory / f"{level}.fna.gz", "rb") as handle:
            result[level] = {
                name: sequence.decode()
                for name, sequence in fragmentation.read_fasta(handle)
            }
    return result


class LevelTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmp = Path(tmp.name)
        self.chr1 = (
            "ACGT" * 5 + "N" * 10 + "GGCC" * 3 + "N" * 5 + "TTAA" + "n" * 12 + "CA"
        )
        self.records = {
            "OZ1.1": self.chr1,
            "OZ2.1": "N" * 12 + "AC" * 10,
            "CAX1.1": "GATTACA",
            "CAX2.1": "CCCC" + "N" * 10 + "GGGG",
            "OZ3.1": "ATATATAT",
        }
        self.report = [
            seq_row("OZ1.1", "assembled-molecule", "1"),
            seq_row("OZ2.1", "assembled-molecule", "2"),
            seq_row("CAX1.1", "unlocalized-scaffold", "1"),
            seq_row("CAX2.1", "unplaced-scaffold", "Un"),
            seq_row(
                "OZ3.1", "assembled-molecule", "MIT", "non-nuclear", "Mitochondrion"
            ),
        ]

    def write_genome(self, records):
        path = self.tmp / "genome.fna.gz"
        with gzip.open(path, "wt") as handle:
            for name, sequence in records.items():
                handle.write(f">{name} description\n")
                for i in range(0, len(sequence), 7):
                    handle.write(sequence[i : i + 7] + "\n")
        return path

    def test_levels_split_gaps_and_exclude_organelles(self):
        stats = fragmentation.build_levels(
            self.write_genome(self.records), self.report, self.tmp / "levels"
        )
        levels = read_levels(self.tmp / "levels")
        self.assertEqual(
            levels["chromosomes"], {"chr1": self.chr1, "chr2": self.records["OZ2.1"]}
        )
        self.assertEqual(
            list(levels["scaffolds"]), ["chr1", "chr2", "CAX1.1", "CAX2.1"]
        )
        # Only N runs of at least MIN_GAP split; a leading gap leaves one piece
        # that keeps its scaffold name.
        self.assertEqual(
            levels["contigs"],
            {
                "chr1_ctg1": "ACGT" * 5,
                "chr1_ctg2": "GGCC" * 3 + "N" * 5 + "TTAA",
                "chr1_ctg3": "CA",
                "chr2": "AC" * 10,
                "CAX1.1": "GATTACA",
                "CAX2.1_ctg1": "CCCC",
                "CAX2.1_ctg2": "GGGG",
            },
        )
        self.assertEqual(stats["organelles_excluded"], 1)
        self.assertEqual(
            stats["levels"]["chromosomes"],
            {
                "sequences": 2,
                "total_length": len(self.chr1) + 32,
                "n50": len(self.chr1),
                "n_bases": 27 + 12,
            },
        )
        self.assertEqual(stats["levels"]["contigs"]["sequences"], 7)
        self.assertEqual(stats["levels"]["contigs"]["n_bases"], 5)
        self.assertEqual(stats["levels"]["contigs"]["n50"], 20)

    def test_fasta_and_report_must_describe_the_same_sequences(self):
        extra = dict(self.records, UNKNOWN="ACGT")
        with self.assertRaisesRegex(ValueError, "absent from the sequence report"):
            fragmentation.build_levels(
                self.write_genome(extra), self.report, self.tmp / "a"
            )
        missing = {k: v for k, v in self.records.items() if k != "CAX1.1"}
        with self.assertRaisesRegex(ValueError, "not in the FASTA"):
            fragmentation.build_levels(
                self.write_genome(missing), self.report, self.tmp / "b"
            )

    def test_duplicate_chromosome_names_fall_back_to_accessions(self):
        names = fragmentation.classify_sequences(
            [
                seq_row("OZ1.1", "assembled-molecule", "1"),
                seq_row("OZ9.1", "assembled-molecule", "1"),
            ]
        )
        self.assertEqual(names["OZ1.1"], ("chromosome", "chr1"))
        self.assertEqual(names["OZ9.1"], ("chromosome", "OZ9.1"))

    def test_unknown_nuclear_role_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unsupported sequence role"):
            fragmentation.classify_sequences([seq_row("X1.1", "alt-scaffold")])


def write_taxdump(directory):
    nodes = [
        (1, 1, "no rank"),
        (10, 1, "phylum"),
        (20, 1, "phylum"),
        (30, 1, "phylum"),
        (100, 10, "genus"),
        (101, 100, "species"),
        (102, 100, "species"),
        (103, 100, "species"),
        (201, 20, "species"),
        (301, 30, "species"),
    ]
    names = {10: "Mollusca", 20: "Chordata", 30: "Rotifera"}
    (directory / "nodes.dmp").write_text(
        "".join(f"{t}\t|\t{p}\t|\t{r}\t|\n" for t, p, r in nodes)
    )
    (directory / "names.dmp").write_text(
        "".join(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for t, n in names.items())
    )
    (directory / "merged.dmp").write_text("999\t|\t102\t|\n")


def report(accession, taxid, **info):
    assembly_info = {
        "assembly_name": f"asm{accession}",
        "assembly_level": "Chromosome",
        "assembly_status": "current",
        "assembly_type": "haploid",
        "release_date": "2025-01-01",
        **info,
    }
    return {
        "accession": accession,
        "organism": {"organism_name": f"Species {taxid}", "tax_id": taxid},
        "assembly_info": assembly_info,
        "assembly_stats": {"number_of_contigs": 10},
    }


class SamplingTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        write_taxdump(Path(tmp.name))
        self.taxonomy = fragmentation.Taxonomy(tmp.name)
        self.reports = [
            report("GCA_1.1", 101),
            report("GCA_2.1", 101, release_date="2026-01-01"),
            report("GCF_2.1", 101, release_date="2026-01-01"),
            report("GCA_3.1", 102, assembly_type="alternate-pseudohaplotype"),
            report("GCA_4.1", 103),
            report("GCA_5.1", 201),
            report("GCA_6.1", 301),
            report("GCA_7.1", 555),
            report("GCA_8.1", 999, assembly_level="Contig"),
            report("GCA_9.1", 102, assembly_status="suppressed"),
        ]
        self.sketches = {"Mollusca", "Chordata"}

    def test_exclusions_and_species_representatives(self):
        rows = fragmentation.select_candidates(
            self.reports, self.taxonomy, self.sketches
        )
        self.assertEqual(
            {row["accession"]: row["status"] for row in rows},
            {
                "GCA_1.1": "excluded:same_species_duplicate",
                "GCA_2.1": "eligible",
                "GCF_2.1": "excluded:not_genbank_accession",
                "GCA_3.1": "excluded:alternate_haplotype",
                "GCA_4.1": "eligible",
                "GCA_5.1": "eligible",
                "GCA_6.1": "excluded:no_mash_sketch",
                "GCA_7.1": "excluded:taxid_not_in_taxdump",
                "GCA_8.1": "excluded:not_chromosome_level",
                "GCA_9.1": "excluded:not_current",
            },
        )

    def test_one_genome_per_phylum_independent_of_input_order(self):
        def draw(reports, seed):
            rows = fragmentation.select_candidates(
                reports, self.taxonomy, self.sketches
            )
            return [
                r["accession"] for r in fragmentation.sample_one_per_phylum(rows, seed)
            ]

        sampled = draw(self.reports, 7)
        self.assertEqual(sampled[0], "GCA_5.1")  # Chordata, the only candidate
        self.assertIn(sampled[1], {"GCA_2.1", "GCA_4.1"})  # Mollusca
        self.assertEqual(len(sampled), 2)
        self.assertEqual(draw(list(reversed(self.reports)), 7), sampled)
        draws = {draw(self.reports, seed)[1] for seed in range(20)}
        self.assertEqual(draws, {"GCA_2.1", "GCA_4.1"})


class MetricTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmp = Path(tmp.name)

    def test_significant_anchor_fraction(self):
        paf = self.tmp / "aln_busco.paf"
        pairs = [("c1", "r1")] * 3 + [("c2", "r1")] + [("c3", "r2")] * 2
        paf.write_text(
            "".join(
                f"{q}\t100\t1\t2\t+\t{t}\t100\t1\t2\t1000\t1000\t60\n" for q, t in pairs
            )
        )
        significance = self.tmp / "significance_results.txt"
        significance.write_text(
            "query\ttarget\traw_p\tcorrected_p\tsignificant\n"
            "c1\tr1\t1e-9\t1e-8\ttrue\n"
            "c2\tr1\t1e-1\t1e0\tfalse\n"
            "c3\tr2\t1e-6\t1e-5\ttrue\n"
        )
        self.assertEqual(
            fragmentation.anchor_metrics(paf, significance),
            {
                "busco_anchors": 6,
                "anchor_sequences": 3,
                "significant_pairs": 2,
                "significant_anchor_fraction": 0.8333,
            },
        )

    def test_busco_summary_tolerates_extra_fields(self):
        directory = self.tmp / "busco_assembly"
        directory.mkdir()
        (
            directory / "short_summary.specific.mollusca_odb12.busco_assembly.txt"
        ).write_text("\tC:95.2%[S:94.0%,D:1.2%],F:1.0%,M:3.8%,n:5295,E:0.5%\n")
        self.assertEqual(
            fragmentation.busco_summary(directory),
            {"C": 95.2, "S": 94.0, "D": 1.2, "F": 1.0, "M": 3.8, "n": 5295.0},
        )

    def test_benchmark_step_groups(self):
        header = "s\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time\n"
        for stem, seconds in (
            ("busco_assembly", "100"),
            ("busco_reference_GCA_1.1", "90"),
            ("busco_to_paf_GCA_1.1", "2"),
            ("aln_GCA_1.1", "50"),
            ("select_references", "1"),
            ("get_reference_GCA_1.1", "30"),
            ("jcvi_synteny", "NA"),
        ):
            (self.tmp / f"{stem}.txt").write_text(
                header + f"{seconds}\t0:00:00\t10\t-\t-\t-\t-\t-\t-\t5\n"
            )
        groups = {
            row["rule"]: (row["group"], row["s"])
            for row in fragmentation.read_benchmarks(self.tmp)
        }
        self.assertEqual(
            groups,
            {
                "busco_assembly": ("busco_assembly", 100.0),
                "busco_reference_GCA_1.1": ("busco_reference", 90.0),
                "busco_to_paf_GCA_1.1": ("other", 2.0),
                "aln_GCA_1.1": ("minimap2", 50.0),
                "select_references": ("mash", 1.0),
                "get_reference_GCA_1.1": ("download", 30.0),
                "jcvi_synteny": ("other", None),
            },
        )

    def test_hobrac_command_rejects_design_changing_options(self):
        row = {"organism_name": "Felimare picta (x)", "taxid": "1"}
        command = fragmentation.hobrac_command(
            row, "in.fna.gz", "out", "20260629", "apptainer", ["--metaeuk"]
        )
        self.assertEqual(command[command.index("-n") + 1], "Felimare picta _x_")
        self.assertEqual(command[-2:], ["--use-apptainer", "--metaeuk"])
        for option in ("--skip-genomic", "-r", "--mash-db-version=latest"):
            with self.assertRaises(ValueError):
                fragmentation.hobrac_command(row, "a", "o", "v", "none", [option])


if __name__ == "__main__":
    unittest.main()
