#!/usr/bin/env python3
"""
Generates multi-species karyotype plots with statistical ALG (Ancestral Linkage
Group) detection using Fisher's exact test with Bonferroni correction.
"""

import argparse
import glob
import os
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple

from scipy.stats import fisher_exact


@dataclass
class BuscoGene:
    busco_id: str
    chromosome: str
    start: int
    end: int


@dataclass
class ALGAssociation:
    chr1: str
    chr2: str
    p_value: float
    color: str
    gene_count: int


# 37-color palette for ALG visualization
ALG_PALETTE = [
    "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
    "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
    "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000",
    "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080",
    "#ff6f61", "#6b5b95", "#88b04b", "#f7cac9", "#92a8d1",
    "#955251", "#b565a7", "#009b77", "#dd4124", "#d65076",
    "#45b8ac", "#efc050", "#5b5ea6", "#9b2335", "#dfcfbe",
    "#55b4b0", "#e15d44"
]


def read_busco_tsv(file_path: str) -> Dict[str, BuscoGene]:
    """
    Parse BUSCO full_table.tsv and return only Complete single-copy genes.

    Args:
        file_path: Path to BUSCO full_table.tsv file

    Returns:
        Dictionary mapping BUSCO ID to BuscoGene
    """
    busco_data = {}
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 5 and parts[1] == "Complete":
                busco_id = parts[0]
                busco_data[busco_id] = BuscoGene(
                    busco_id=busco_id,
                    chromosome=parts[2],
                    start=int(parts[3]),
                    end=int(parts[4])
                )
    return busco_data


def detect_algs_pairwise(
    busco1: Dict[str, BuscoGene],
    busco2: Dict[str, BuscoGene],
    alpha: float = 0.01,
    min_genes: int = 5
) -> Tuple[List[ALGAssociation], Dict[str, str]]:
    """
    Detect ALG associations using Fisher's exact test with Bonferroni correction.

    Args:
        busco1: BUSCO data for species 1
        busco2: BUSCO data for species 2
        alpha: Significance level before correction
        min_genes: Minimum genes in a chr pair to test

    Returns:
        Tuple of (list of significant ALG associations, dict mapping gene to color)
    """
    # Find common BUSCO IDs
    common_ids = set(busco1.keys()) & set(busco2.keys())
    if not common_ids:
        return [], {}

    # Count genes per chromosome pair
    chr_pair_counts = defaultdict(int)
    chr1_counts = defaultdict(int)
    chr2_counts = defaultdict(int)

    for busco_id in common_ids:
        chr1 = busco1[busco_id].chromosome
        chr2 = busco2[busco_id].chromosome
        chr_pair_counts[(chr1, chr2)] += 1
        chr1_counts[chr1] += 1
        chr2_counts[chr2] += 1

    total_genes = len(common_ids)

    # Collect all chromosome pairs that meet minimum gene threshold
    testable_pairs = [
        (chr1, chr2, count)
        for (chr1, chr2), count in chr_pair_counts.items()
        if count >= min_genes
    ]

    if not testable_pairs:
        return [], {}

    # Perform Fisher's exact test for each pair
    results = []
    for chr1, chr2, observed in testable_pairs:
        # Build 2x2 contingency table
        # [observed, chr1_other], [chr2_other, neither]
        chr1_other = chr1_counts[chr1] - observed
        chr2_other = chr2_counts[chr2] - observed
        neither = total_genes - chr1_counts[chr1] - chr2_counts[chr2] + observed

        table = [[observed, chr1_other], [chr2_other, neither]]
        _, p_value = fisher_exact(table, alternative='greater')
        results.append((chr1, chr2, p_value, observed))

    # Apply Bonferroni correction
    num_tests = len(results)
    p_threshold = alpha / num_tests

    # Filter significant associations and assign colors
    significant = []
    color_idx = 0
    for chr1, chr2, p_value, gene_count in results:
        if p_value < p_threshold:
            color = ALG_PALETTE[color_idx % len(ALG_PALETTE)]
            significant.append(ALGAssociation(
                chr1=chr1, chr2=chr2, p_value=p_value,
                color=color, gene_count=gene_count
            ))
            color_idx += 1

    # Build gene-to-color mapping
    gene_colors = {}
    for busco_id in common_ids:
        chr1 = busco1[busco_id].chromosome
        chr2 = busco2[busco_id].chromosome
        color = "lightgrey"
        for alg in significant:
            if alg.chr1 == chr1 and alg.chr2 == chr2:
                color = alg.color
                break
        gene_colors[busco_id] = color

    return significant, gene_colors


def generate_bed_file(
    busco_data: Dict[str, BuscoGene],
    species_name: str,
    output_path: str
) -> None:
    """
    Generate a BED file for JCVI from BUSCO data.

    Args:
        busco_data: Dictionary of BUSCO genes
        species_name: Species name prefix for gene IDs
        output_path: Output BED file path
    """
    # Sort by chromosome and position
    genes = sorted(
        busco_data.values(),
        key=lambda g: (g.chromosome, g.start)
    )

    with open(output_path, 'w') as f:
        for gene in genes:
            # BED format: chr, start, end, gene_id, score, strand
            gene_id = f"{species_name}_{gene.busco_id}"
            f.write(f"{gene.chromosome}\t{gene.start}\t{gene.end}\t{gene_id}\t0\t+\n")


def generate_links_file(
    busco1: Dict[str, BuscoGene],
    busco2: Dict[str, BuscoGene],
    gene_colors: Dict[str, str],
    sp1: str,
    sp2: str,
    output_path: str
) -> None:
    """
    Generate JCVI links file between two species.

    Args:
        busco1: BUSCO data for species 1
        busco2: BUSCO data for species 2
        gene_colors: Dictionary mapping BUSCO ID to color
        sp1: Species 1 name
        sp2: Species 2 name
        output_path: Output links file path
    """
    common_ids = set(busco1.keys()) & set(busco2.keys())

    with open(output_path, 'w') as f:
        for busco_id in sorted(common_ids):
            gene1_id = f"{sp1}_{busco_id}"
            gene2_id = f"{sp2}_{busco_id}"
            color = gene_colors.get(busco_id, "lightgrey")
            f.write(f"{gene1_id}\t{gene2_id}\t{color}\n")


def get_chromosome_order(
    busco_data: Dict[str, BuscoGene]
) -> List[str]:
    """
    Get chromosomes sorted by total gene span (largest first).

    Args:
        busco_data: BUSCO gene data

    Returns:
        List of chromosome names sorted by size
    """
    chr_spans = defaultdict(lambda: [float('inf'), 0])
    for gene in busco_data.values():
        chr_spans[gene.chromosome][0] = min(chr_spans[gene.chromosome][0], gene.start)
        chr_spans[gene.chromosome][1] = max(chr_spans[gene.chromosome][1], gene.end)

    # Sort by span (largest first)
    sorted_chrs = sorted(
        chr_spans.keys(),
        key=lambda c: chr_spans[c][1] - chr_spans[c][0],
        reverse=True
    )
    return sorted_chrs


def generate_seqids_file(
    species_data: List[Tuple[str, Dict[str, BuscoGene]]],
    output_path: str
) -> None:
    """
    Generate JCVI seqids file with chromosome order for each species.

    Args:
        species_data: List of (species_name, busco_data) tuples
        output_path: Output seqids file path
    """
    with open(output_path, 'w') as f:
        for species_name, busco_data in species_data:
            chromosomes = get_chromosome_order(busco_data)
            f.write(",".join(chromosomes) + "\n")


def generate_layouts_file(
    species_beds: List[Tuple[str, str]],
    links_files: List[str],
    output_path: str
) -> None:
    """
    Generate JCVI layouts file.

    Args:
        species_beds: List of (species_name, bed_file_path) tuples
        links_files: List of links file paths
        output_path: Output layouts file path
    """
    num_species = len(species_beds)

    with open(output_path, 'w') as f:
        f.write("# y, xstart, xend, rotation, color, label, va, bed\n")
        f.write(f"#{'-' * 60}\n")

        # Calculate y positions (evenly spaced from 0.9 to 0.1)
        y_positions = [0.9 - (i * 0.8 / max(num_species - 1, 1))
                       for i in range(num_species)]

        for i, (species_name, bed_path) in enumerate(species_beds):
            bed_basename = os.path.basename(bed_path)
            y = y_positions[i]
            va = "top" if i == 0 else "bottom"
            f.write(f"{y:.2f}, 0.1, 0.9, 0, black, {species_name}, {va}, {bed_basename}\n")

        f.write("\n# edges\n")
        for links_file in links_files:
            links_basename = os.path.basename(links_file)
            f.write(f"e, 0, 1, {links_basename}\n")


def save_alg_associations(
    algs: List[ALGAssociation],
    sp1: str,
    sp2: str,
    output_path: str
) -> None:
    """
    Save ALG associations to a TSV file for inspection.

    Args:
        algs: List of ALG associations
        sp1: Species 1 name
        sp2: Species 2 name
        output_path: Output TSV file path
    """
    with open(output_path, 'a') as f:
        if os.path.getsize(output_path) == 0:
            f.write("species1\tspecies2\tchr1\tchr2\tp_value\tgene_count\tcolor\n")
        for alg in algs:
            f.write(f"{sp1}\t{sp2}\t{alg.chr1}\t{alg.chr2}\t"
                    f"{alg.p_value:.2e}\t{alg.gene_count}\t{alg.color}\n")


def get_species_order(
    assembly_name: str,
    assembly_busco_path: str,
    accession_order_file: str,
    manual_refs: str,
    busco_refs_pattern: str
) -> List[Tuple[str, str]]:
    """
    Determine species order and collect BUSCO paths.

    Args:
        assembly_name: Name for the assembly
        assembly_busco_path: Path to assembly BUSCO full_table.tsv
        accession_order_file: Path to MASH-ordered accessions file
        manual_refs: Semicolon-separated manual reference paths (or empty)
        busco_refs_pattern: Pattern to find reference BUSCO directories

    Returns:
        List of (species_name, busco_path) tuples in display order
    """
    species = [(assembly_name, assembly_busco_path)]

    if manual_refs:
        # Manual mode: use order given, label = base filename
        for ref_path in manual_refs.split(";"):
            name = os.path.splitext(os.path.basename(ref_path))[0]
            busco_path = find_busco_table(busco_refs_pattern, name)
            if busco_path:
                species.append((name, busco_path))
    else:
        # MASH-discovered: order from file (already sorted by distance)
        with open(accession_order_file) as f:
            for line in f:
                accession = line.strip()
                if accession:
                    busco_path = find_busco_table(busco_refs_pattern, accession)
                    if busco_path:
                        species.append((accession, busco_path))

    return species


def find_busco_table(pattern: str, accession: str) -> str:
    """
    Find the BUSCO full_table.tsv for a given accession.

    Args:
        pattern: Glob pattern containing {accession} placeholder
        accession: Accession ID to search for

    Returns:
        Path to full_table.tsv or empty string if not found
    """
    search_pattern = pattern.replace("{accession}", accession)
    matches = glob.glob(search_pattern)
    if matches:
        return matches[0]
    return ""


def run(
    assembly_busco: str,
    busco_refs: List[str],
    accession_order: str,
    manual_refs: str,
    output_dir: str,
    assembly_name: str = "assembly"
) -> Dict[str, str]:
    """
    Main entry point for JCVI synteny analysis.

    Args:
        assembly_busco: Path to assembly BUSCO full_table.tsv
        busco_refs: List of reference BUSCO directory paths
        accession_order: Path to MASH-ordered accessions file
        manual_refs: Semicolon-separated manual reference paths
        output_dir: Output directory path
        assembly_name: Name for the assembly in plots

    Returns:
        Dictionary with paths to generated files
    """
    os.makedirs(output_dir, exist_ok=True)

    # Build mapping from accession to BUSCO path
    ref_busco_paths = {}
    for ref_dir in busco_refs:
        # Extract accession from directory name (busco_reference_ACCESSION)
        dirname = os.path.basename(ref_dir.rstrip('/'))
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

    # Load all BUSCO data
    species_busco = []
    for name, path in species_order:
        busco_data = read_busco_tsv(path)
        species_busco.append((name, busco_data))

    # Generate BED files
    bed_files = []
    for name, busco_data in species_busco:
        bed_path = os.path.join(output_dir, f"{name}.bed")
        generate_bed_file(busco_data, name, bed_path)
        bed_files.append((name, bed_path))

    # Detect ALGs and generate links for consecutive species pairs
    alg_output = os.path.join(output_dir, "alg_associations.tsv")
    open(alg_output, 'w').close()  # Create empty file

    links_files = []
    for i in range(len(species_busco) - 1):
        sp1_name, sp1_busco = species_busco[i]
        sp2_name, sp2_busco = species_busco[i + 1]

        # Detect ALGs
        algs, gene_colors = detect_algs_pairwise(sp1_busco, sp2_busco)

        # Save ALG associations
        save_alg_associations(algs, sp1_name, sp2_name, alg_output)

        # Generate links file
        links_path = os.path.join(output_dir, f"links.{sp1_name}.{sp2_name}.simple")
        generate_links_file(sp1_busco, sp2_busco, gene_colors,
                            sp1_name, sp2_name, links_path)
        links_files.append(links_path)

    # Generate seqids file
    seqids_path = os.path.join(output_dir, "seqids")
    generate_seqids_file(species_busco, seqids_path)

    # Generate layouts file
    layouts_path = os.path.join(output_dir, "layouts")
    generate_layouts_file(bed_files, links_files, layouts_path)

    return {
        "seqids": seqids_path,
        "layouts": layouts_path,
        "alg_associations": alg_output,
        "bed_files": [p for _, p in bed_files],
        "links_files": links_files
    }


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Generate JCVI karyotype files from BUSCO data with ALG detection"
    )
    parser.add_argument(
        '--busco_assembly',
        required=True,
        help="Path to assembly BUSCO full_table.tsv"
    )
    parser.add_argument(
        '--busco_references',
        nargs='+',
        required=True,
        help="Paths to reference BUSCO directories"
    )
    parser.add_argument(
        '--accession_order',
        required=True,
        help="Path to file with accessions in order (MASH distance)"
    )
    parser.add_argument(
        '--manual_refs',
        default="",
        help="Semicolon-separated manual reference paths"
    )
    parser.add_argument(
        '--output_dir',
        required=True,
        help="Output directory"
    )
    parser.add_argument(
        '--assembly_name',
        default="assembly",
        help="Name for assembly in plots"
    )

    args = parser.parse_args()

    run(
        assembly_busco=args.busco_assembly,
        busco_refs=args.busco_references,
        accession_order=args.accession_order,
        manual_refs=args.manual_refs,
        output_dir=args.output_dir,
        assembly_name=args.assembly_name
    )


if __name__ == "__main__":
    main()
