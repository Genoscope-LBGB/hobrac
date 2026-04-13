import argparse
import glob
import os
from typing import Dict, List, Tuple

from .coloring import apply_custom_colors, apply_custom_colors_with_algs
from .io import parse_custom_colors, read_busco_tsv, read_fasta_sizes
from .models import (
    DEFAULT_COLOR,
    BuscoGene,
    ChromosomeAssociation,
)
from .ordering import (
    get_chromosome_order,
    get_chromosome_order_by_gravity,
    get_chromosome_order_by_span,
)
from .statistics import detect_algs_transitive


def generate_bed_file(
    busco_data: Dict[str, BuscoGene], species_name: str, output_path: str
) -> None:
    """
    Generate a BED file for JCVI from BUSCO data.

    Args:
        busco_data: Dictionary of BUSCO genes
        species_name: Species name prefix for gene IDs
        output_path: Output BED file path
    """
    # Sort by chromosome and position (use min of start/end for sorting)
    genes = sorted(
        busco_data.values(), key=lambda g: (g.chromosome, min(g.start, g.end))
    )

    with open(output_path, "w") as f:
        for gene in genes:
            # BED format: chr, start, end, gene_id, score, strand
            # BED requires start <= end, so swap if necessary
            # Also determine strand based on original coordinate order
            if gene.start <= gene.end:
                start, end, strand = gene.start, gene.end, "+"
            else:
                start, end, strand = gene.end, gene.start, "-"
            gene_id = f"{species_name}_{gene.busco_id}"
            f.write(f"{gene.chromosome}\t{start}\t{end}\t{gene_id}\t0\t{strand}\n")


def generate_links_file(
    busco1: Dict[str, BuscoGene],
    busco2: Dict[str, BuscoGene],
    gene_colors: Dict[str, str],
    sp1: str,
    sp2: str,
    output_path: str,
    hide_non_significant: bool = False,
) -> None:
    """
    Generate a JCVI .simple file with one line per BUSCO gene pair.

    Each line represents a single ortholog and uses its ``gene_colors`` entry
    as the link color. Because ``gene_colors`` is derived from the single,
    run-wide ``gene_to_chain`` mapping, this guarantees that a gene keeps the
    same color across every .simple file of a single run — no block-level
    majority vote can override its assignment.

    JCVI karyotype expects the 6-column format:
    ``startGene1 endGene1 startGene2 endGene2 score orientation``

    Args:
        busco1: BUSCO data for species 1
        busco2: BUSCO data for species 2
        gene_colors: Dictionary mapping BUSCO ID to color
        sp1: Species 1 name
        sp2: Species 2 name
        output_path: Output .simple file path
        hide_non_significant: If True, skip genes colored lightgrey
    """
    common_ids = set(busco1.keys()) & set(busco2.keys())
    if not common_ids:
        open(output_path, "w").close()
        return

    # Non-significant (grey) genes are written first so that JCVI draws the
    # significant ribbons on top of them. Ties are broken by gene id for
    # deterministic output.
    sorted_ids = sorted(
        common_ids,
        key=lambda gid: (
            gene_colors.get(gid, DEFAULT_COLOR) != DEFAULT_COLOR,
            gid,
        ),
    )

    with open(output_path, "w") as f:
        for busco_id in sorted_ids:
            color = gene_colors.get(busco_id, DEFAULT_COLOR)
            if hide_non_significant and color == DEFAULT_COLOR:
                continue
            gene1 = f"{sp1}_{busco_id}"
            gene2 = f"{sp2}_{busco_id}"
            f.write(f"{color}*{gene1}\t{gene1}\t{gene2}\t{gene2}\t1\t+\n")


def generate_seqids_file(
    species_data: List[Tuple[str, Dict[str, BuscoGene]]],
    output_path: str,
    use_gravity_ordering: bool = False,
    assembly_fasta_sizes: Dict[str, int] = None,
) -> None:
    """
    Generate JCVI seqids file with chromosome order for each species.

    Args:
        species_data: List of (species_name, busco_data) tuples
        output_path: Output seqids file path
        use_gravity_ordering: If True, order subsequent species by gravity
                              with the species above them
        assembly_fasta_sizes: Dictionary mapping assembly sequence names to sizes.
                              Used for ordering the first species (assembly).
    """
    with open(output_path, "w") as f:
        previous_order: List[str] = []
        previous_busco: Dict[str, BuscoGene] = {}

        for i, (species_name, busco_data) in enumerate(species_data):
            if i == 0:
                # First species (assembly): use fasta size ordering
                chromosomes = get_chromosome_order(assembly_fasta_sizes, busco_data)
            elif use_gravity_ordering:
                # Use gravity ordering with the species above
                chromosomes = get_chromosome_order_by_gravity(
                    previous_order, busco_data, previous_busco
                )
            else:
                # Gravity disabled: use span ordering for references
                chromosomes = get_chromosome_order_by_span(busco_data)

            f.write(",".join(chromosomes) + "\n")

            # Store for next iteration
            previous_order = chromosomes
            previous_busco = busco_data


def generate_layouts_file(
    species_beds: List[Tuple[str, str]], links_files: List[str], output_path: str
) -> None:
    """
    Generate JCVI layouts file.

    Args:
        species_beds: List of (species_name, bed_file_path) tuples
        links_files: List of links file paths
        output_path: Output layouts file path
    """
    num_species = len(species_beds)

    with open(output_path, "w") as f:
        f.write("# y, xstart, xend, rotation, color, label, va, bed\n")
        f.write(f"#{'-' * 60}\n")

        # Calculate y positions (evenly spaced from 0.9 to 0.1)
        y_positions = [
            0.9 - (i * 0.8 / max(num_species - 1, 1)) for i in range(num_species)
        ]

        for i, (species_name, bed_path) in enumerate(species_beds):
            bed_basename = os.path.basename(bed_path)
            y = y_positions[i]
            va = "top" if i == 0 else "bottom"
            f.write(
                f"{y:.2f},\t0.1,\t0.9,\t0,\tblack,\t{species_name},\t{va},\t{bed_basename}\n"
            )

        f.write("\n# edges\n")
        for i, links_file in enumerate(links_files):
            links_basename = os.path.basename(links_file)
            f.write(f"e, {i}, {i + 1}, {links_basename}\n")


def save_chromosome_associations(
    associations: List[ChromosomeAssociation], output_path: str
) -> None:
    """
    Save chromosome associations to a TSV file with header.

    Args:
        associations: List of chromosome associations
        output_path: Output TSV file path
    """
    with open(output_path, "w") as f:
        f.write(
            "species1\tspecies2\tchr1\tchr2\tp_value\t"
            "corrected_p_value\tgene_count\tsignificant\n"
        )
        for a in associations:
            significant = "accepted" if a.significant else "rejected"
            f.write(
                f"{a.species1}\t{a.species2}\t{a.chr1}\t{a.chr2}\t"
                f"{a.p_value:.2e}\t{a.corrected_p_value:.2e}\t"
                f"{a.gene_count}\t{significant}\n"
            )


def run(
    assembly_busco: str,
    assembly_fasta: str,
    busco_refs: List[str],
    accession_order: str,
    manual_refs: str,
    output_dir: str,
    assembly_name: str = "assembly",
    use_gravity_ordering: bool = True,
    min_busco_genes: int = 0,
    custom_color_file: str = "",
    custom_names: str = "",
    hide_non_significant: bool = False,
    skip_alg: bool = False,
) -> Dict[str, str]:
    """
    Main entry point for JCVI synteny analysis.

    Args:
        assembly_busco: Path to assembly BUSCO full_table.tsv
        assembly_fasta: Path to the assembly fasta file
        busco_refs: List of reference BUSCO directory paths
        accession_order: Path to MASH-ordered accessions file
        manual_refs: Semicolon-separated manual reference paths
        output_dir: Output directory path
        assembly_name: Name for the assembly in plots
        use_gravity_ordering: If True, order chromosomes by gravity for
                              diagonal alignment patterns
        min_busco_genes: Minimum complete BUSCO genes required per sequence
        custom_color_file: Path to custom color file for gene coloring.
        custom_names: Comma-separated custom names for JCVI tracks. If one name
                      is provided, it applies only to the assembly. If multiple
                      names are provided, the count must match species count.
        hide_non_significant: If True, hide links between chromosome pairs
                              without significant associations.
        skip_alg: If True, skip ALG statistical testing when using custom
                  colors. Only effective with custom_color_file.

    Returns:
        Dictionary with paths to generated files
    """
    os.makedirs(output_dir, exist_ok=True)

    # Read assembly fasta sizes for chromosome ordering
    assembly_fasta_sizes = read_fasta_sizes(assembly_fasta)

    # Load custom colors if provided
    custom_colors = {}
    if custom_color_file:
        custom_colors = parse_custom_colors(custom_color_file)

    # Build mapping from accession to BUSCO path
    ref_busco_paths = {}
    for ref_dir in busco_refs:
        # Extract accession from directory name (busco_reference_ACCESSION)
        dirname = os.path.basename(ref_dir.rstrip("/"))
        if dirname.startswith("busco_reference_"):
            accession = dirname.replace("busco_reference_", "")
        else:
            accession = dirname

        # Find full_table.tsv within the directory
        table_pattern = os.path.join(ref_dir, "run*", "full_table.tsv")
        matches = glob.glob(table_pattern)
        if matches:
            ref_busco_paths[accession] = matches[0]

    # Determine species order
    species_order = []
    species_order.append((assembly_name, assembly_busco))

    if manual_refs:
        for ref_path in manual_refs.split(";"):
            if ref_path.strip():
                name = os.path.splitext(os.path.basename(ref_path))[0]
                if name in ref_busco_paths:
                    species_order.append((name, ref_busco_paths[name]))
    else:
        with open(accession_order) as f:
            for line in f:
                accession = line.strip()
                if accession and accession in ref_busco_paths:
                    species_order.append((accession, ref_busco_paths[accession]))

    # Apply custom names if provided
    if custom_names:
        names_list = [n.strip() for n in custom_names.split(",")]
        if len(names_list) == 1:
            # Single name: apply only to assembly
            if species_order:
                species_order[0] = (names_list[0], species_order[0][1])
        else:
            # Multiple names: must match species count
            if len(names_list) != len(species_order):
                raise ValueError(
                    f"Number of custom names ({len(names_list)}) must equal "
                    f"number of species ({len(species_order)}): 1 assembly + "
                    f"{len(species_order) - 1} references"
                )
            species_order = [
                (names_list[i], species_order[i][1]) for i in range(len(species_order))
            ]

    # Load all BUSCO data
    species_busco = []
    for name, path in species_order:
        busco_data = read_busco_tsv(path, min_busco_genes)
        species_busco.append((name, busco_data))

    # Generate BED files
    bed_files = []
    for name, busco_data in species_busco:
        bed_path = os.path.join(output_dir, f"{name}.bed")
        generate_bed_file(busco_data, name, bed_path)
        bed_files.append((name, bed_path))

    # Detect ALGs and generate links for consecutive species pairs
    associations_output = os.path.join(output_dir, "chromosome_associations.tsv")

    # Phase 1: Transitive ALG detection across all species
    # Skip only when custom colors are provided AND skip_alg is True
    all_chromosome_associations: List[ChromosomeAssociation] = []
    gene_to_chain: Dict[str, int] = {}
    chain_colors: Dict[int, str] = {}
    if not (custom_colors and skip_alg):
        _, _, gene_to_chain, chain_colors, all_chromosome_associations = (
            detect_algs_transitive(species_busco)
        )
    save_chromosome_associations(all_chromosome_associations, associations_output)

    # Phase 2: Generate outputs for each pair
    links_files = []
    for i in range(len(species_busco) - 1):
        sp1_name, sp1_busco = species_busco[i]
        sp2_name, sp2_busco = species_busco[i + 1]

        if custom_colors and skip_alg:
            gene_colors = apply_custom_colors(sp1_busco, sp2_busco, custom_colors)
        else:
            gene_colors = apply_custom_colors_with_algs(
                sp1_busco,
                sp2_busco,
                gene_to_chain,
                chain_colors,
                custom_colors,
            )

        links_path = os.path.join(output_dir, f"links.{sp1_name}.{sp2_name}.simple")
        generate_links_file(
            sp1_busco,
            sp2_busco,
            gene_colors,
            sp1_name,
            sp2_name,
            links_path,
            hide_non_significant=hide_non_significant,
        )
        links_files.append(links_path)

    seqids_path = os.path.join(output_dir, "seqids")
    generate_seqids_file(
        species_busco, seqids_path, use_gravity_ordering, assembly_fasta_sizes
    )

    layouts_path = os.path.join(output_dir, "layouts")
    generate_layouts_file(bed_files, links_files, layouts_path)

    return {
        "seqids": seqids_path,
        "layouts": layouts_path,
        "chromosome_associations": associations_output,
        "bed_files": [p for _, p in bed_files],
        "links_files": links_files,
    }


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Generate JCVI karyotype files from BUSCO data with ALG detection"
    )
    parser.add_argument(
        "--busco_assembly", required=True, help="Path to assembly BUSCO full_table.tsv"
    )
    parser.add_argument(
        "--assembly-fasta", required=True, help="Path to the assembly fasta file"
    )
    parser.add_argument(
        "--busco_references",
        nargs="+",
        required=True,
        help="Paths to reference BUSCO directories",
    )
    parser.add_argument(
        "--accession_order",
        required=True,
        help="Path to file with accessions in order (MASH distance)",
    )
    parser.add_argument(
        "--manual_refs", default="", help="Semicolon-separated manual reference paths"
    )
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument(
        "--assembly_name", default="assembly", help="Name for assembly in plots"
    )
    parser.add_argument(
        "--min-busco-genes",
        type=int,
        default=0,
        help="Minimum complete BUSCO genes required per sequence (default: 0)",
    )
    parser.add_argument(
        "--jcvi-custom-colors",
        default="",
        help="Path to custom color file (tab-separated: BUSCO_ID, R,G,B, ALG_NAME). "
        "By default, ALG statistical testing still runs to determine significance; "
        "use --skip-alg to disable it. Genes not in the file will be shown in grey.",
    )
    parser.add_argument(
        "--jcvi-names",
        default="",
        help="Comma-separated custom names for JCVI tracks. If one name is provided, "
        "it applies only to the assembly. If multiple names are provided, the count "
        "must equal 1 (assembly) + number of references, in order.",
    )
    parser.add_argument(
        "--hide-non-significant",
        action="store_true",
        help="Hide links between chromosome pairs without significant associations. "
        "This produces a cleaner plot showing only ALG-related synteny.",
    )
    parser.add_argument(
        "--skip-alg",
        action="store_true",
        help="Skip ALG statistical testing when using custom colors. All genes in the "
        "color file get their custom color; unlisted genes are grey. Only effective "
        "with --jcvi-custom-colors.",
    )

    args = parser.parse_args()

    run(
        assembly_busco=args.busco_assembly,
        assembly_fasta=args.assembly_fasta,
        busco_refs=args.busco_references,
        accession_order=args.accession_order,
        manual_refs=args.manual_refs,
        output_dir=args.output_dir,
        assembly_name=args.assembly_name,
        min_busco_genes=args.min_busco_genes,
        custom_color_file=args.jcvi_custom_colors,
        custom_names=args.jcvi_names,
        hide_non_significant=args.hide_non_significant,
        skip_alg=args.skip_alg,
    )


if __name__ == "__main__":
    main()
