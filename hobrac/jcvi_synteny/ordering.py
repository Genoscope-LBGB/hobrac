import math
from collections import defaultdict
from typing import Dict, List, Tuple

from .models import BuscoGene


def get_chromosome_order(
    fasta_sizes: Dict[str, int], busco_data: Dict[str, BuscoGene]
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
        busco_chromosomes, key=lambda c: fasta_sizes.get(c, 0), reverse=True
    )
    return sorted_chrs


def get_chromosome_order_by_span(busco_data: Dict[str, BuscoGene]) -> List[str]:
    """
    Get chromosomes sorted by total gene span (largest first).

    Args:
        busco_data: BUSCO gene data

    Returns:
        List of chromosome names sorted by gene span
    """
    chr_spans = defaultdict(lambda: [float("inf"), 0])
    for gene in busco_data.values():
        chr_spans[gene.chromosome][0] = min(chr_spans[gene.chromosome][0], gene.start)
        chr_spans[gene.chromosome][1] = max(chr_spans[gene.chromosome][1], gene.end)

    sorted_chrs = sorted(
        chr_spans.keys(), key=lambda c: chr_spans[c][1] - chr_spans[c][0], reverse=True
    )
    return sorted_chrs


def calculate_gravity_scores(
    target_busco: Dict[str, BuscoGene], query_busco: Dict[str, BuscoGene]
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
        distance = math.sqrt(target_span**2 + query_span**2)

        # Gravity formula: (1 + distance)^2
        gravity = (1 + distance) ** 2

        chr_pair = (target_gene.chromosome, query_gene.chromosome)
        gravity_scores[chr_pair] += gravity

    return dict(gravity_scores)


def get_chromosome_order_by_gravity(
    target_order: List[str],
    query_busco: Dict[str, BuscoGene],
    target_busco: Dict[str, BuscoGene],
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
            group.sort(key=lambda q: query_positions.get(q, float("inf")))
            ordered_queries.extend(group)
            matched_queries.update(group)

    # Add unmatched chromosomes at the end, sorted by span
    all_query_chrs = set(gene.chromosome for gene in query_busco.values())
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
