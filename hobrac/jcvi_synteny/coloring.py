from typing import Dict, List, Tuple

from .models import DEFAULT_COLOR, BuscoGene


def apply_custom_colors(
    busco1: Dict[str, BuscoGene],
    busco2: Dict[str, BuscoGene],
    custom_colors: Dict[str, str],
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
        gene_colors[busco_id] = custom_colors.get(busco_id, DEFAULT_COLOR)
    return gene_colors


def _chain_covers_pair(
    chain: List[Tuple[str, str]], sp1_name: str, sp2_name: str
) -> bool:
    """Return True if *chain* contains both species of the pair."""
    species_in_chain = {node[0] for node in chain}
    return sp1_name in species_in_chain and sp2_name in species_in_chain


def apply_custom_colors_with_algs(
    species1_busco: Dict[str, BuscoGene],
    species2_busco: Dict[str, BuscoGene],
    gene_to_chain: Dict[str, int],
    chain_colors: Dict[int, str],
    custom_colors: Dict[str, str],
    chains: List[List[Tuple[str, str]]],
    sp1_name: str,
    sp2_name: str,
) -> Dict[str, str]:
    """
    Determine link colors using the gene's selected display chain.

    Only the selected chain is considered, not the union of all compatible
    chains. A link is colored only when that chain contains both species;
    direct pairwise significance is not re-tested here. Otherwise the link
    is lightgrey, even when the gene has a custom color.

    Genes covered by their display chain get their custom color when listed,
    the chain palette color when no custom mapping is provided, or lightgrey
    when a custom mapping is provided but does not list the gene.

    Args:
        species1_busco: BUSCO data for species 1
        species2_busco: BUSCO data for species 2
        gene_to_chain: BUSCO gene ID to selected display chain ID (-1 = none)
        chain_colors: Dict mapping chain ID to hex color
        custom_colors: Dictionary mapping BUSCO ID to hex color
        chains: Retained chains from ALG detection
        sp1_name: Name of species 1 in the current pair
        sp2_name: Name of species 2 in the current pair

    Returns:
        Dictionary mapping BUSCO ID to color
    """
    chains_covering_pair = {
        i
        for i, chain in enumerate(chains)
        if _chain_covers_pair(chain, sp1_name, sp2_name)
    }

    common_ids = set(species1_busco.keys()) & set(species2_busco.keys())
    gene_colors = {}

    for busco_id in common_ids:
        chain_id = gene_to_chain.get(busco_id, -1)
        if chain_id < 0 or chain_id not in chains_covering_pair:
            gene_colors[busco_id] = DEFAULT_COLOR
        elif custom_colors:
            gene_colors[busco_id] = custom_colors.get(busco_id, DEFAULT_COLOR)
        else:
            gene_colors[busco_id] = chain_colors[chain_id]

    return gene_colors
