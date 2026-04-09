from collections import defaultdict
from typing import Dict, List, Set, Tuple

from .models import BuscoGene, PairwiseAssociation


def enumerate_chains(
    pairwise_associations: List[PairwiseAssociation],
) -> List[List[Tuple[str, str]]]:
    """
    Enumerate all maximal chromosome chains through the DAG of significant
    pairwise associations across consecutive species.

    Each chain is a maximal path where edges go from species i to species i+1.

    Args:
        pairwise_associations: List of significant PairwiseAssociation objects

    Returns:
        List of chains sorted deterministically, where each chain is a list
        of (species, chromosome) tuples in species order.
    """
    if not pairwise_associations:
        return []

    outgoing: Dict[Tuple[str, str], List[Tuple[str, str]]] = defaultdict(list)
    has_incoming: Set[Tuple[str, str]] = set()

    for assoc in pairwise_associations:
        src = (assoc.species1, assoc.chr1)
        dst = (assoc.species2, assoc.chr2)
        outgoing[src].append(dst)
        has_incoming.add(dst)

    for successors in outgoing.values():
        successors.sort()

    roots = sorted(set(outgoing) - has_incoming)

    chains: List[List[Tuple[str, str]]] = []

    def dfs(node: Tuple[str, str], path: List[Tuple[str, str]]) -> None:
        successors = outgoing.get(node)
        if not successors:
            chains.append(list(path))
            return
        for succ in successors:
            path.append(succ)
            dfs(succ, path)
            path.pop()

    for root in roots:
        dfs(root, [root])

    chains.sort()
    return chains


def build_gene_chain_mapping(
    species_busco: List[Tuple[str, Dict[str, BuscoGene]]],
    chains: List[List[Tuple[str, str]]],
) -> Dict[str, int]:
    """
    Assign each BUSCO gene to a chain based on its chromosome path across
    all species.

    A gene matches a chain if, at every species covered by the chain where the
    gene is present, the gene's chromosome equals the chain's chromosome.

    Args:
        species_busco: Ordered list of (species_name, busco_data) tuples
        chains: Chains from enumerate_chains

    Returns:
        Dict mapping BUSCO gene ID to chain index (-1 if no matching chain).
        Only genes present in at least two consecutive species are included.
    """
    if not species_busco:
        return {}

    species_index = {name: i for i, (name, _) in enumerate(species_busco)}

    chain_paths = [
        {species_index[sp]: chrom for sp, chrom in chain} for chain in chains
    ]

    eligible = set()
    for i in range(len(species_busco) - 1):
        eligible |= species_busco[i][1].keys() & species_busco[i + 1][1].keys()

    gene_mapping: Dict[str, int] = {}
    chain_counts: Dict[int, int] = defaultdict(int)
    ambiguous: List[Tuple[str, List[int]]] = []

    for gene_id in sorted(eligible):
        gene_chroms: Dict[int, str] = {}
        for pos, (_, busco_data) in enumerate(species_busco):
            if gene_id in busco_data:
                gene_chroms[pos] = busco_data[gene_id].chromosome

        matches = [
            cid
            for cid, cpath in enumerate(chain_paths)
            if _gene_matches_chain(gene_chroms, cpath)
        ]

        if len(matches) == 1:
            gene_mapping[gene_id] = matches[0]
            chain_counts[matches[0]] += 1
        elif len(matches) > 1:
            ambiguous.append((gene_id, matches))
        else:
            gene_mapping[gene_id] = -1

    for gene_id, matches in ambiguous:
        best = min(matches, key=lambda c: (-chain_counts[c], c))
        gene_mapping[gene_id] = best
        chain_counts[best] += 1

    return gene_mapping


def _gene_matches_chain(
    gene_chroms: Dict[int, str], chain_path: Dict[int, str]
) -> bool:
    """Return True if the gene is compatible with the chain at all shared positions."""
    compared = False
    for pos, chrom in chain_path.items():
        if pos in gene_chroms:
            if gene_chroms[pos] != chrom:
                return False
            compared = True
    return compared
