from typing import Dict

from .models import BuscoGene, DEFAULT_COLOR


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


def apply_custom_colors_with_algs(
    species1_busco: Dict[str, BuscoGene],
    species2_busco: Dict[str, BuscoGene],
    gene_to_chain: Dict[str, int],
    chain_colors: Dict[int, str],
    custom_colors: Dict[str, str],
) -> Dict[str, str]:
    """
    Determine link colors for a species pair, gated by chain significance.

    Non-significant genes (not on any chain) are always lightgrey, regardless
    of whether they appear in custom_colors. Significant genes (on a chain)
    get their custom color when listed in custom_colors, the chain palette
    color when no custom file is provided, or lightgrey when a custom file is
    provided but the gene is not listed.

    Args:
        species1_busco: BUSCO data for species 1
        species2_busco: BUSCO data for species 2
        gene_to_chain: Dict mapping BUSCO gene ID to chain ID (-1 = no chain)
        chain_colors: Dict mapping chain ID to hex color
        custom_colors: Dictionary mapping BUSCO ID to hex color

    Returns:
        Dictionary mapping BUSCO ID to color
    """
    common_ids = set(species1_busco.keys()) & set(species2_busco.keys())
    gene_colors = {}

    for busco_id in common_ids:
        chain_id = gene_to_chain.get(busco_id, -1)
        if chain_id < 0:
            gene_colors[busco_id] = DEFAULT_COLOR
        elif custom_colors:
            gene_colors[busco_id] = custom_colors.get(busco_id, DEFAULT_COLOR)
        else:
            gene_colors[busco_id] = chain_colors[chain_id]

    return gene_colors
