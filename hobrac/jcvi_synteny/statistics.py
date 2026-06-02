from collections import defaultdict
from typing import Dict, List, Tuple

from scipy.stats import fisher_exact

from .chains import build_gene_chain_mapping, enumerate_chains, validate_chains
from .models import (
    ALG_PALETTE,
    BuscoGene,
    ChromosomeAssociation,
    PairwiseAssociation,
)


def detect_algs_pairwise_raw(
    busco1: Dict[str, BuscoGene],
    busco2: Dict[str, BuscoGene],
    species1: str,
    species2: str,
    alpha: float = 0.01,
    min_genes: int = 5,
) -> Tuple[List[PairwiseAssociation], List[ChromosomeAssociation]]:
    """
    Run Fisher's exact test, return significant associations and all tested pairs.

    Args:
        busco1: BUSCO data for species 1
        busco2: BUSCO data for species 2
        species1: Name of species 1
        species2: Name of species 2
        alpha: Significance level before correction
        min_genes: Minimum genes in a chr pair to test

    Returns:
        Tuple of:
        - List of significant PairwiseAssociation objects
        - List of all tested ChromosomeAssociation objects (significant and not)
    """
    common_ids = set(busco1.keys()) & set(busco2.keys())
    if not common_ids:
        return [], []

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
        return [], []

    results = []
    for chr1, chr2, observed in testable_pairs:
        chr1_other = chr1_counts[chr1] - observed
        chr2_other = chr2_counts[chr2] - observed
        neither = total_genes - chr1_counts[chr1] - chr2_counts[chr2] + observed

        table = [[observed, chr1_other], [chr2_other, neither]]
        _, p_value = fisher_exact(table, alternative="greater")
        results.append((chr1, chr2, p_value, observed))

    num_tests = len(results)

    significant = []
    all_tested = []
    for chr1, chr2, p_value, gene_count in results:
        corrected_p_value = min(p_value * num_tests, 1.0)
        is_significant = corrected_p_value <= alpha

        all_tested.append(
            ChromosomeAssociation(
                species1=species1,
                species2=species2,
                chr1=chr1,
                chr2=chr2,
                p_value=p_value,
                corrected_p_value=corrected_p_value,
                gene_count=gene_count,
                significant=is_significant,
            )
        )

        if is_significant:
            significant.append(
                PairwiseAssociation(
                    species1=species1,
                    species2=species2,
                    chr1=chr1,
                    chr2=chr2,
                    p_value=p_value,
                    gene_count=gene_count,
                )
            )

    return significant, all_tested


def detect_algs_transitive(
    species_busco: List[Tuple[str, Dict[str, BuscoGene]]],
    alpha: float = 0.01,
    min_genes: int = 5,
    min_chain_genes: int = 5,
    permissive_alg: bool = False,
) -> Tuple[
    List[PairwiseAssociation],
    List[List[Tuple[str, str]]],
    Dict[str, int],
    Dict[int, str],
    List[ChromosomeAssociation],
]:
    """
    Detect ALGs with consistent colors across all species using chain-based
    grouping.

    Args:
        species_busco: List of (species_name, busco_data) tuples
        alpha: Significance level before correction
        min_genes: Minimum genes in a chr pair to test
        min_chain_genes: Minimum BUSCO genes a chain must be supported by to
            survive. Chains with fewer genes are dropped before sub-chain
            pruning, so sub-chains contained in a dropped chain can
            re-emerge. Default 5.

    Returns:
        Tuple of:
        - all_associations: List of PairwiseAssociation objects
        - chains: List of chromosome chains
        - gene_to_chain: Dict mapping BUSCO gene ID to chain ID
        - chain_colors: Dict mapping chain ID to hex color
        - all_chromosome_associations: List of all tested ChromosomeAssociation objects
    """
    all_associations: List[PairwiseAssociation] = []
    all_chromosome_associations: List[ChromosomeAssociation] = []
    for i in range(len(species_busco)):
        for j in range(i + 1, len(species_busco)):
            sp1_name, sp1_busco = species_busco[i]
            sp2_name, sp2_busco = species_busco[j]
            significant, tested = detect_algs_pairwise_raw(
                sp1_busco, sp2_busco, sp1_name, sp2_name, alpha, min_genes
            )
            all_associations.extend(significant)
            all_chromosome_associations.extend(tested)

    chains = enumerate_chains(
        all_associations, species_busco, min_chain_genes=min_chain_genes
    )
    if all_chromosome_associations:
        chains = validate_chains(
            chains,
            all_chromosome_associations,
            permissive=permissive_alg,
            min_chain_genes=min_chain_genes,
            species_busco=species_busco,
        )

    gene_to_chain = build_gene_chain_mapping(species_busco, chains)
    chain_colors = {i: ALG_PALETTE[i % len(ALG_PALETTE)] for i in range(len(chains))}

    return (
        all_associations,
        chains,
        gene_to_chain,
        chain_colors,
        all_chromosome_associations,
    )
