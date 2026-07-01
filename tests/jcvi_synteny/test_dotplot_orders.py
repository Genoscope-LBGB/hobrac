"""Tests for the ALG-dotplot axis-order files (generate_seqids_file +
save_dotplot_axis_orders)."""

from hobrac.jcvi_synteny.models import BuscoGene
from hobrac.jcvi_synteny.output import (
    generate_seqids_file,
    save_dotplot_axis_orders,
)


def _gene(chromosome, start=0, end=100):
    return BuscoGene(busco_id="", chromosome=chromosome, start=start, end=end)


def test_generate_seqids_returns_orders_matching_file(tmp_path):
    # Assembly ordered by fasta size (chr1 largest), one reference.
    species_data = [
        ("Felimare_picta", {"g1": _gene("chr1"), "g2": _gene("chr2")}),
        ("GCA_123", {"g1": _gene("scfA"), "g2": _gene("scfB")}),
    ]
    assembly_sizes = {"chr1": 5000, "chr2": 1000}
    seqids = tmp_path / "seqids"

    orders = generate_seqids_file(
        species_data, str(seqids), use_gravity_ordering=True,
        assembly_fasta_sizes=assembly_sizes,
    )

    # One order list per species, and the assembly is size-ordered.
    assert len(orders) == 2
    assert orders[0] == ["chr1", "chr2"]

    # The returned orders match exactly what was written to the seqids file.
    file_lines = [line.strip().split(",") for line in seqids.read_text().splitlines()]
    assert file_lines == orders


def test_save_dotplot_axis_orders_keys_by_accession(tmp_path):
    species_keys = ["Felimare_picta", "GCA_123"]
    orders = [["chr1", "chr2"], ["scfA", "scfB"]]

    save_dotplot_axis_orders(species_keys, orders, str(tmp_path))

    orders_dir = tmp_path / "dotplot_orders"
    assert (orders_dir / "Felimare_picta.order").read_text().strip() == "chr1,chr2"
    assert (orders_dir / "GCA_123.order").read_text().strip() == "scfA,scfB"
