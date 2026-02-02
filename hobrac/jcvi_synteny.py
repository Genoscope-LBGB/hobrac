#!/usr/bin/env python3
"""
Generates multi-species karyotype plots with statistical ALG (Ancestral Linkage
Group) detection using Fisher's exact test with Bonferroni correction.
"""

import argparse
import glob
import os
from collections import defaultdict, deque
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple

from scipy.stats import fisher_exact


@dataclass
class BuscoGene:
    busco_id: str
    chromosome: str
    start: int
    end: int


@dataclass
class PairwiseAssociation:
    """Raw pairwise significant association before ALG grouping."""
    species1: str
    species2: str
    chr1: str
    chr2: str
    p_value: float
    gene_count: int


@dataclass
class ALGAssociation:
    chr1: str
    chr2: str
    p_value: float
    color: str
    gene_count: int
    alg_id: int = -1


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


def read_fasta_sizes(fasta_path: str) -> Dict[str, int]:
    """
    Read sequence sizes from a fasta file.

    Args:
        fasta_path: Path to the fasta file

    Returns:
        Dictionary mapping sequence names to their lengths
    """
    sizes = {}
    current_name = None
    current_length = 0

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name is not None:
                    sizes[current_name] = current_length
                current_name = line[1:].split()[0]
                current_length = 0
            else:
                current_length += len(line)

        if current_name is not None:
            sizes[current_name] = current_length

    return sizes


def parse_custom_colors(color_file: str) -> Dict[str, str]:
    """
    Parse custom color file and convert RGB values to hex colors.

    Color file format (tab or space separated):
    BUSCO_ID     R,G,B        ALG_NAME
    10000at6447  141,78,106   A2

    Args:
        color_file: Path to the custom color file

    Returns:
        Dictionary mapping BUSCO ID to hex color string
    """
    custom_colors = {}
    with open(color_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            busco_id = parts[0]
            rgb_str = parts[1]
            try:
                r, g, b = map(int, rgb_str.split(','))
                hex_color = f"#{r:02x}{g:02x}{b:02x}"
                custom_colors[busco_id] = hex_color
            except (ValueError, IndexError):
                # Skip malformed lines
                continue
    return custom_colors


def apply_custom_colors(
    busco1: Dict[str, BuscoGene],
    busco2: Dict[str, BuscoGene],
    custom_colors: Dict[str, str]
) -> Dict[str, str]:
    """
    Apply custom colors to genes common between two species.

    Args:
        busco1: BUSCO data for species 1
        busco2: BUSCO data for species 2
        custom_colors: Dictionary mapping BUSCO ID to hex color

    Returns:
        Dictionary mapping BUSCO ID to color (custom or lightgrey)
    """
    common_ids = set(busco1.keys()) & set(busco2.keys())
    gene_colors = {}
    for busco_id in common_ids:
        gene_colors[busco_id] = custom_colors.get(busco_id, "lightgrey")
    return gene_colors


def build_alg_graph(
    pairwise_associations: List[PairwiseAssociation]
) -> Dict[Tuple[str, str], Set[Tuple[str, str]]]:
    """
    Build adjacency list from significant associations.

    Args:
        pairwise_associations: List of PairwiseAssociation objects

    Returns:
        Graph where nodes are (species, chromosome) tuples
    """
    graph: Dict[Tuple[str, str], Set[Tuple[str, str]]] = defaultdict(set)
    for assoc in pairwise_associations:
        node1 = (assoc.species1, assoc.chr1)
        node2 = (assoc.species2, assoc.chr2)
        graph[node1].add(node2)
        graph[node2].add(node1)
    return dict(graph)


def find_connected_components(
    graph: Dict[Tuple[str, str], Set[Tuple[str, str]]]
) -> List[Set[Tuple[str, str]]]:
    """
    Find connected components

    Args:
        graph: Adjacency list representation of the graph

    Returns:
        List of component sets, where each set contains (species, chr) tuples
    """
    if not graph:
        return []

    visited: Set[Tuple[str, str]] = set()
    components: List[Set[Tuple[str, str]]] = []

    for start_node in graph:
        if start_node in visited:
            continue

        component: Set[Tuple[str, str]] = set()
        queue = deque([start_node])

        while queue:
            node = queue.popleft()
            if node in visited:
                continue
            visited.add(node)
            component.add(node)

            for neighbor in graph.get(node, set()):
                if neighbor not in visited:
                    queue.append(neighbor)

        components.append(component)

    return components


def read_busco_tsv(
    file_path: str,
    min_busco_genes: int = 0
) -> Dict[str, BuscoGene]:
    """
    Parse BUSCO full_table.tsv and return only Complete single-copy genes.

    Args:
        file_path: Path to BUSCO full_table.tsv file
        min_busco_genes: Minimum complete BUSCO genes required per sequence.
                         Sequences with fewer genes are excluded.

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

    if min_busco_genes > 0:
        busco_data = filter_by_min_genes(busco_data, min_busco_genes)

    return busco_data


def filter_by_min_genes(
    busco_data: Dict[str, BuscoGene],
    min_genes: int
) -> Dict[str, BuscoGene]:
    """
    Filter BUSCO data to keep only sequences with at least min_genes.

    Args:
        busco_data: Dictionary mapping BUSCO ID to BuscoGene
        min_genes: Minimum number of genes required per sequence

    Returns:
        Filtered dictionary with only genes from qualifying sequences
    """
    seq_counts: Dict[str, int] = defaultdict(int)
    for gene in busco_data.values():
        seq_counts[gene.chromosome] += 1

    valid_seqs = {seq for seq, count in seq_counts.items() if count >= min_genes}

    return {
        busco_id: gene
        for busco_id, gene in busco_data.items()
        if gene.chromosome in valid_seqs
    }


def detect_algs_pairwise_raw(
    busco1: Dict[str, BuscoGene],
    busco2: Dict[str, BuscoGene],
    species1: str,
    species2: str,
    alpha: float = 0.01,
    min_genes: int = 5
) -> List[PairwiseAssociation]:
    """
    Run Fisher's exact test, return raw associations without colors.

    Args:
        busco1: BUSCO data for species 1
        busco2: BUSCO data for species 2
        species1: Name of species 1
        species2: Name of species 2
        alpha: Significance level before correction
        min_genes: Minimum genes in a chr pair to test

    Returns:
        List of significant PairwiseAssociation objects
    """
    common_ids = set(busco1.keys()) & set(busco2.keys())
    if not common_ids:
        return []

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

    testable_pairs = [
        (chr1, chr2, count)
        for (chr1, chr2), count in chr_pair_counts.items()
        if count >= min_genes
    ]

    if not testable_pairs:
        return []

    results = []
    for chr1, chr2, observed in testable_pairs:
        chr1_other = chr1_counts[chr1] - observed
        chr2_other = chr2_counts[chr2] - observed
        neither = total_genes - chr1_counts[chr1] - chr2_counts[chr2] + observed

        table = [[observed, chr1_other], [chr2_other, neither]]
        _, p_value = fisher_exact(table, alternative='greater')
        results.append((chr1, chr2, p_value, observed))

    num_tests = len(results)
    p_threshold = alpha / num_tests

    significant = []
    for chr1, chr2, p_value, gene_count in results:
        if p_value < p_threshold:
            significant.append(PairwiseAssociation(
                species1=species1,
                species2=species2,
                chr1=chr1,
                chr2=chr2,
                p_value=p_value,
                gene_count=gene_count
            ))

    return significant


def detect_algs_transitive(
    species_busco: List[Tuple[str, Dict[str, BuscoGene]]],
    alpha: float = 0.01,
    min_genes: int = 5
) -> Tuple[List[PairwiseAssociation], Dict[Tuple[str, str], int], Dict[int, str]]:
    """
    Detect ALGs with consistent colors across all species.

    Args:
        species_busco: List of (species_name, busco_data) tuples
        alpha: Significance level before correction
        min_genes: Minimum genes in a chr pair to test

    Returns:
        Tuple of:
        - all_associations: List of PairwiseAssociation objects
        - chr_to_alg: Dict mapping (species, chr) to alg_id
        - alg_colors: Dict mapping alg_id to color
    """
    all_associations: List[PairwiseAssociation] = []
    for i in range(len(species_busco) - 1):
        sp1_name, sp1_busco = species_busco[i]
        sp2_name, sp2_busco = species_busco[i + 1]
        associations = detect_algs_pairwise_raw(
            sp1_busco, sp2_busco, sp1_name, sp2_name, alpha, min_genes
        )
        all_associations.extend(associations)

    graph = build_alg_graph(all_associations)
    components = find_connected_components(graph)

    chr_to_alg: Dict[Tuple[str, str], int] = {}
    alg_colors: Dict[int, str] = {}

    for alg_id, component in enumerate(components):
        color = ALG_PALETTE[alg_id % len(ALG_PALETTE)]
        alg_colors[alg_id] = color
        for node in component:
            chr_to_alg[node] = alg_id

    return all_associations, chr_to_alg, alg_colors


def build_gene_colors_from_algs(
    species1_busco: Dict[str, BuscoGene],
    species2_busco: Dict[str, BuscoGene],
    species1: str,
    species2: str,
    chr_to_alg: Dict[Tuple[str, str], int],
    alg_colors: Dict[int, str],
    significant_associations: List[PairwiseAssociation]
) -> Dict[str, str]:
    """
    Build gene-to-color mapping for a species pair using ALG membership.

    Only genes on chromosome pairs with significant associations are colored.
    Genes on non-significant chromosome pairs are colored grey.

    Args:
        species1_busco: BUSCO data for species 1
        species2_busco: BUSCO data for species 2
        species1: Name of species 1
        species2: Name of species 2
        chr_to_alg: Dict mapping (species, chr) to alg_id
        alg_colors: Dict mapping alg_id to color
        significant_associations: List of significant associations for this pair

    Returns:
        Dictionary mapping BUSCO ID to color
    """
    # Build set of significant chromosome pairs
    significant_pairs = {
        (assoc.chr1, assoc.chr2) for assoc in significant_associations
    }

    common_ids = set(species1_busco.keys()) & set(species2_busco.keys())
    gene_colors = {}

    for busco_id in common_ids:
        chr1 = species1_busco[busco_id].chromosome
        chr2 = species2_busco[busco_id].chromosome

        # Only color if this specific chromosome pair has a significant association
        if (chr1, chr2) in significant_pairs:
            node1 = (species1, chr1)
            alg_id = chr_to_alg.get(node1)
            if alg_id is not None:
                gene_colors[busco_id] = alg_colors[alg_id]
            else:
                gene_colors[busco_id] = "lightgrey"
        else:
            gene_colors[busco_id] = "lightgrey"

    return gene_colors


def detect_algs_pairwise(
    busco1: Dict[str, BuscoGene],
    busco2: Dict[str, BuscoGene],
    alpha: float = 0.01,
    min_genes: int = 5
) -> Tuple[List[ALGAssociation], Dict[str, str]]:
    """
    Detect ALG associations using Fisher's exact test with Bonferroni correction.

    This is a wrapper for backward compatibility (single pair use case).

    Args:
        busco1: BUSCO data for species 1
        busco2: BUSCO data for species 2
        alpha: Significance level before correction
        min_genes: Minimum genes in a chr pair to test

    Returns:
        Tuple of (list of significant ALG associations, dict mapping gene to color)
    """
    raw_associations = detect_algs_pairwise_raw(
        busco1, busco2, "sp1", "sp2", alpha, min_genes
    )

    significant = []
    for i, assoc in enumerate(raw_associations):
        color = ALG_PALETTE[i % len(ALG_PALETTE)]
        significant.append(ALGAssociation(
            chr1=assoc.chr1,
            chr2=assoc.chr2,
            p_value=assoc.p_value,
            color=color,
            gene_count=assoc.gene_count,
            alg_id=i
        ))

    common_ids = set(busco1.keys()) & set(busco2.keys())
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
    # Sort by chromosome and position (use min of start/end for sorting)
    genes = sorted(
        busco_data.values(),
        key=lambda g: (g.chromosome, min(g.start, g.end))
    )

    with open(output_path, 'w') as f:
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
    max_gap: int = 10,
    hide_non_significant: bool = False
) -> None:
    """
    Generate JCVI .simple file with synteny blocks between two species.

    JCVI karyotype expects 6-column format:
    startGene1 endGene1 startGene2 endGene2 score orientation

    Args:
        busco1: BUSCO data for species 1
        busco2: BUSCO data for species 2
        gene_colors: Dictionary mapping BUSCO ID to color
        sp1: Species 1 name
        sp2: Species 2 name
        output_path: Output .simple file path
        max_gap: Maximum rank gap to consider genes as part of same block
        hide_non_significant: If True, skip blocks with no significant association
    """
    common_ids = set(busco1.keys()) & set(busco2.keys())
    if not common_ids:
        # Write empty file
        open(output_path, 'w').close()
        return

    # Build chromosome-sorted gene lists for each species
    # and assign positional ranks within each chromosome
    def get_chr_ranks(busco_data: Dict[str, BuscoGene]) -> Dict[str, Tuple[str, int]]:
        """Return dict: busco_id -> (chromosome, rank_on_chromosome)"""
        # Group by chromosome
        chr_genes = defaultdict(list)
        for busco_id, gene in busco_data.items():
            if busco_id in common_ids:
                chr_genes[gene.chromosome].append((busco_id, min(gene.start, gene.end)))

        # Sort each chromosome by position and assign ranks
        result = {}
        for chrom, genes in chr_genes.items():
            genes.sort(key=lambda x: x[1])
            for rank, (busco_id, _) in enumerate(genes):
                result[busco_id] = (chrom, rank)
        return result

    ranks1 = get_chr_ranks(busco1)
    ranks2 = get_chr_ranks(busco2)

    # Group orthologs by chromosome pair
    chr_pair_orthologs = defaultdict(list)
    for busco_id in common_ids:
        if busco_id in ranks1 and busco_id in ranks2:
            chr1, rank1 = ranks1[busco_id]
            chr2, rank2 = ranks2[busco_id]
            chr_pair_orthologs[(chr1, chr2)].append((busco_id, rank1, rank2))

    # For each chromosome pair, identify syntenic blocks
    blocks = []
    for (chr1, chr2), orthologs in chr_pair_orthologs.items():
        if len(orthologs) < 2:
            # Single gene - still output as a block
            if orthologs:
                busco_id, rank1, rank2 = orthologs[0]
                blocks.append({
                    'start1': busco_id, 'end1': busco_id,
                    'start2': busco_id, 'end2': busco_id,
                    'score': 1, 'orientation': '+',
                    'color': gene_colors.get(busco_id, 'lightgrey')
                })
            continue

        # Sort by position in species 1
        orthologs.sort(key=lambda x: x[1])

        # Detect blocks: consecutive in sp1, check direction in sp2
        current_block = [orthologs[0]]
        for i in range(1, len(orthologs)):
            busco_id, rank1, rank2 = orthologs[i]
            prev_busco, prev_rank1, prev_rank2 = current_block[-1]

            # Check if this gene continues the block
            gap1 = rank1 - prev_rank1
            gap2 = abs(rank2 - prev_rank2)

            if gap1 <= max_gap and gap2 <= max_gap:
                # Continue current block
                current_block.append(orthologs[i])
            else:
                # Save current block and start new one
                if current_block:
                    blocks.append(_make_block(current_block, gene_colors))
                current_block = [orthologs[i]]

        # Don't forget last block
        if current_block:
            blocks.append(_make_block(current_block, gene_colors))

    # Write .simple file in JCVI format
    # Color prefix format: "color*gene_name" (e.g., "#ff0000*gene1")
    with open(output_path, 'w') as f:
        for block in blocks:
            color = block.get('color', 'lightgrey')
            # Skip non-significant blocks if requested
            if hide_non_significant and color == 'lightgrey':
                continue
            # Add color prefix if not default grey
            color_prefix = f"{color}*" if color != 'lightgrey' else ""
            gene1_start = f"{color_prefix}{sp1}_{block['start1']}"
            gene1_end = f"{sp1}_{block['end1']}"
            gene2_start = f"{sp2}_{block['start2']}"
            gene2_end = f"{sp2}_{block['end2']}"
            f.write(f"{gene1_start}\t{gene1_end}\t{gene2_start}\t{gene2_end}\t"
                    f"{block['score']}\t{block['orientation']}\n")


def _make_block(
    orthologs: List[Tuple[str, int, int]],
    gene_colors: Dict[str, str]
) -> Dict:
    """
    Create a synteny block from a list of orthologs.

    Args:
        orthologs: List of (busco_id, rank1, rank2) tuples, sorted by rank1
        gene_colors: Dictionary mapping BUSCO ID to color

    Returns:
        Block dictionary with start/end genes, score, orientation, color
    """
    first = orthologs[0]
    last = orthologs[-1]

    # Determine orientation based on rank2 direction
    rank2_first = first[2]
    rank2_last = last[2]
    orientation = '+' if rank2_last >= rank2_first else '-'

    # For inverted blocks, swap start2/end2 so they're in correct order
    if orientation == '+':
        start2, end2 = first[0], last[0]
    else:
        start2, end2 = last[0], first[0]

    # Get most common color in block (for ALG coloring)
    color_counts = defaultdict(int)
    for busco_id, _, _ in orthologs:
        color = gene_colors.get(busco_id, 'lightgrey')
        color_counts[color] += 1
    dominant_color = max(color_counts.keys(), key=lambda c: color_counts[c])

    return {
        'start1': first[0],
        'end1': last[0],
        'start2': start2,
        'end2': end2,
        'score': len(orthologs),
        'orientation': orientation,
        'color': dominant_color
    }


def get_chromosome_order(
    fasta_sizes: Dict[str, int],
    busco_data: Dict[str, BuscoGene]
) -> List[str]:
    """
    Get chromosomes sorted by fasta sequence size (largest first).

    Only chromosomes that have BUSCO genes are included in the output.

    Args:
        fasta_sizes: Dictionary mapping sequence names to their lengths
        busco_data: BUSCO gene data (used to filter to BUSCO-containing sequences)

    Returns:
        List of chromosome names sorted by fasta size (largest first)
    """
    busco_chromosomes = {gene.chromosome for gene in busco_data.values()}

    sorted_chrs = sorted(
        busco_chromosomes,
        key=lambda c: fasta_sizes.get(c, 0),
        reverse=True
    )
    return sorted_chrs


def get_chromosome_order_by_span(
    busco_data: Dict[str, BuscoGene]
) -> List[str]:
    """
    Get chromosomes sorted by total gene span (largest first).

    Args:
        busco_data: BUSCO gene data

    Returns:
        List of chromosome names sorted by gene span
    """
    chr_spans = defaultdict(lambda: [float('inf'), 0])
    for gene in busco_data.values():
        chr_spans[gene.chromosome][0] = min(chr_spans[gene.chromosome][0], gene.start)
        chr_spans[gene.chromosome][1] = max(chr_spans[gene.chromosome][1], gene.end)

    sorted_chrs = sorted(
        chr_spans.keys(),
        key=lambda c: chr_spans[c][1] - chr_spans[c][0],
        reverse=True
    )
    return sorted_chrs


def calculate_gravity_scores(
    target_busco: Dict[str, BuscoGene],
    query_busco: Dict[str, BuscoGene]
) -> Dict[Tuple[str, str], float]:
    """
    Calculate gravity scores between chromosome pairs using shared BUSCO genes.

    Gravity score = sum of (1 + distance)^2 for all shared genes,
    where distance is the Euclidean distance based on gene spans.

    Args:
        target_busco: BUSCO data for target species (species above in plot)
        query_busco: BUSCO data for query species (species to order)

    Returns:
        Dictionary mapping (target_chr, query_chr) to gravity score
    """
    import math

    common_ids = set(target_busco.keys()) & set(query_busco.keys())
    if not common_ids:
        return {}

    gravity_scores: Dict[Tuple[str, str], float] = defaultdict(float)

    for busco_id in common_ids:
        target_gene = target_busco[busco_id]
        query_gene = query_busco[busco_id]

        # Calculate distance using gene spans
        target_span = abs(target_gene.end - target_gene.start)
        query_span = abs(query_gene.end - query_gene.start)
        distance = math.sqrt(target_span ** 2 + query_span ** 2)

        # Gravity formula: (1 + distance)^2
        gravity = (1 + distance) ** 2

        chr_pair = (target_gene.chromosome, query_gene.chromosome)
        gravity_scores[chr_pair] += gravity

    return dict(gravity_scores)


def get_chromosome_order_by_gravity(
    target_order: List[str],
    query_busco: Dict[str, BuscoGene],
    target_busco: Dict[str, BuscoGene]
) -> List[str]:
    """
    Order query chromosomes by gravity with target chromosomes.

    Each query chromosome is matched to the target chromosome with highest
    cumulative gravity. Within each target group, chromosomes are ordered
    by their earliest alignment position on the target.

    Args:
        target_order: Ordered list of target chromosome names
        query_busco: BUSCO data for query species
        target_busco: BUSCO data for target species

    Returns:
        Ordered list of query chromosome names
    """
    gravity_scores = calculate_gravity_scores(target_busco, query_busco)

    if not gravity_scores:
        # Fall back to span ordering if no shared genes
        return get_chromosome_order_by_span(query_busco)

    # Find best matching target for each query chromosome
    query_to_best_target: Dict[str, str] = {}
    query_to_gravity: Dict[str, float] = {}

    # Group scores by query chromosome
    query_chr_scores: Dict[str, Dict[str, float]] = defaultdict(dict)
    for (target_chr, query_chr), score in gravity_scores.items():
        query_chr_scores[query_chr][target_chr] = score

    # Find best target for each query
    for query_chr, target_scores in query_chr_scores.items():
        best_target = max(target_scores.keys(), key=lambda t: target_scores[t])
        query_to_best_target[query_chr] = best_target
        query_to_gravity[query_chr] = target_scores[best_target]

    # Group query chromosomes by their best-matching target
    target_to_queries: Dict[str, List[str]] = defaultdict(list)
    for query_chr, target_chr in query_to_best_target.items():
        target_to_queries[target_chr].append(query_chr)

    # Calculate earliest alignment position for each query on its best target
    query_positions: Dict[str, int] = {}
    common_ids = set(target_busco.keys()) & set(query_busco.keys())

    for busco_id in common_ids:
        target_gene = target_busco[busco_id]
        query_gene = query_busco[busco_id]
        query_chr = query_gene.chromosome

        if query_chr in query_to_best_target:
            if query_to_best_target[query_chr] == target_gene.chromosome:
                pos = min(target_gene.start, target_gene.end)
                if query_chr not in query_positions or pos < query_positions[query_chr]:
                    query_positions[query_chr] = pos

    # Build ordered list following target order
    ordered_queries: List[str] = []
    matched_queries: set = set()

    for target_chr in target_order:
        if target_chr in target_to_queries:
            # Sort queries in this group by position on target
            group = target_to_queries[target_chr]
            group.sort(key=lambda q: query_positions.get(q, float('inf')))
            ordered_queries.extend(group)
            matched_queries.update(group)

    # Add unmatched chromosomes at the end, sorted by span
    all_query_chrs = set(
        gene.chromosome for gene in query_busco.values()
    )
    unmatched = all_query_chrs - matched_queries

    if unmatched:
        # Calculate spans for unmatched chromosomes
        chr_spans: Dict[str, int] = defaultdict(int)
        for gene in query_busco.values():
            if gene.chromosome in unmatched:
                span = abs(gene.end - gene.start)
                chr_spans[gene.chromosome] = max(chr_spans[gene.chromosome], span)

        unmatched_sorted = sorted(unmatched, key=lambda c: chr_spans[c], reverse=True)
        ordered_queries.extend(unmatched_sorted)

    return ordered_queries


def generate_seqids_file(
    species_data: List[Tuple[str, Dict[str, BuscoGene]]],
    output_path: str,
    use_gravity_ordering: bool = False,
    assembly_fasta_sizes: Dict[str, int] = None
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
    with open(output_path, 'w') as f:
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
        for i, links_file in enumerate(links_files):
            links_basename = os.path.basename(links_file)
            f.write(f"e, {i}, {i+1}, {links_basename}\n")


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
            f.write("species1\tspecies2\tchr1\tchr2\tp_value\tgene_count\tcolor\talg_id\n")
        for alg in algs:
            f.write(f"{sp1}\t{sp2}\t{alg.chr1}\t{alg.chr2}\t"
                    f"{alg.p_value:.2e}\t{alg.gene_count}\t{alg.color}\t{alg.alg_id}\n")


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
    hide_non_significant: bool = False
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
        custom_color_file: Path to custom color file. When provided, ALG
                           statistical test is skipped and colors are applied
                           directly from the file.
        custom_names: Comma-separated custom names for JCVI tracks. If one name
                      is provided, it applies only to the assembly. If multiple
                      names are provided, the count must match species count.
        hide_non_significant: If True, hide links between chromosome pairs
                              without significant associations.

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
                (names_list[i], species_order[i][1])
                for i in range(len(species_order))
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
    alg_output = os.path.join(output_dir, "alg_associations.tsv")
    open(alg_output, 'w').close()  # Create empty file

    # Phase 1: Transitive ALG detection across all species (if not using custom colors)
    all_associations: List[PairwiseAssociation] = []
    chr_to_alg: Dict[Tuple[str, str], int] = {}
    alg_colors: Dict[int, str] = {}
    if not custom_colors:
        all_associations, chr_to_alg, alg_colors = detect_algs_transitive(species_busco)

    # Phase 2: Generate outputs for each pair
    links_files = []
    for i in range(len(species_busco) - 1):
        sp1_name, sp1_busco = species_busco[i]
        sp2_name, sp2_busco = species_busco[i + 1]

        if custom_colors:
            # Use custom colors directly, skip statistical analysis
            gene_colors = apply_custom_colors(sp1_busco, sp2_busco, custom_colors)
            algs: List[ALGAssociation] = []
        else:
            # Filter associations for this pair
            pair_associations = [
                a for a in all_associations
                if a.species1 == sp1_name and a.species2 == sp2_name
            ]
            # Build gene colors from ALG membership
            gene_colors = build_gene_colors_from_algs(
                sp1_busco, sp2_busco, sp1_name, sp2_name, chr_to_alg, alg_colors,
                pair_associations
            )
            # Convert to ALGAssociation with alg_id
            algs = [
                ALGAssociation(
                    chr1=a.chr1,
                    chr2=a.chr2,
                    p_value=a.p_value,
                    color=alg_colors.get(chr_to_alg.get((a.species1, a.chr1), -1), "lightgrey"),
                    gene_count=a.gene_count,
                    alg_id=chr_to_alg.get((a.species1, a.chr1), -1)
                )
                for a in pair_associations
            ]

        save_alg_associations(algs, sp1_name, sp2_name, alg_output)

        links_path = os.path.join(output_dir, f"links.{sp1_name}.{sp2_name}.simple")
        generate_links_file(sp1_busco, sp2_busco, gene_colors,
                            sp1_name, sp2_name, links_path,
                            hide_non_significant=hide_non_significant)
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
        '--assembly-fasta',
        required=True,
        help="Path to the assembly fasta file"
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
    parser.add_argument(
        '--min-busco-genes',
        type=int,
        default=0,
        help="Minimum complete BUSCO genes required per sequence (default: 0)"
    )
    parser.add_argument(
        '--jcvi-custom-colors',
        default="",
        help="Path to custom color file (tab-separated: BUSCO_ID, R,G,B, ALG_NAME). "
             "When provided, the ALG statistical test is disabled and colors are "
             "applied directly from the file. Genes not in the file will be shown "
             "in grey."
    )
    parser.add_argument(
        '--jcvi-names',
        default="",
        help="Comma-separated custom names for JCVI tracks. If one name is provided, "
             "it applies only to the assembly. If multiple names are provided, the count "
             "must equal 1 (assembly) + number of references, in order."
    )
    parser.add_argument(
        '--hide-non-significant',
        action='store_true',
        help="Hide links between chromosome pairs without significant associations. "
             "This produces a cleaner plot showing only ALG-related synteny."
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
    )


if __name__ == "__main__":
    main()
