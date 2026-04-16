from collections import defaultdict
from typing import Dict, List, Set, Tuple

from .models import BuscoGene, PairwiseAssociation


def _eligible_genes(
    species_busco: List[Tuple[str, Dict[str, BuscoGene]]],
) -> Set[str]:
    """Return gene IDs present in at least two consecutive species."""
    eligible: Set[str] = set()
    for i in range(len(species_busco) - 1):
        eligible |= species_busco[i][1].keys() & species_busco[i + 1][1].keys()
    return eligible


def enumerate_chains(
    pairwise_associations: List[PairwiseAssociation],
    species_busco: List[Tuple[str, Dict[str, BuscoGene]]],
) -> List[List[Tuple[str, str]]]:
    """
    Enumerate chromosome chains using gene-level evidence to determine which
    paths through significant associations actually exist.

    For each eligible gene, walks its chromosome path across species and
    extracts maximal contiguous sub-paths where all consecutive edges are
    significant associations. Only chains supported by real genes are returned.

    Args:
        pairwise_associations: List of significant PairwiseAssociation objects
        species_busco: Ordered list of (species_name, busco_data) tuples

    Returns:
        List of chains sorted deterministically, where each chain is a list
        of (species, chromosome) tuples in species order.
    """
    if not pairwise_associations or not species_busco:
        return []

    sig_edges: Set[Tuple[Tuple[str, str], Tuple[str, str]]] = set()
    for assoc in pairwise_associations:
        sig_edges.add(((assoc.species1, assoc.chr1), (assoc.species2, assoc.chr2)))

    eligible = _eligible_genes(species_busco)

    observed_chains: Set[Tuple[Tuple[str, str], ...]] = set()

    for gene_id in eligible:
        gene_chroms: Dict[int, Tuple[str, str]] = {}
        for pos, (sp_name, busco_data) in enumerate(species_busco):
            if gene_id in busco_data:
                gene_chroms[pos] = (sp_name, busco_data[gene_id].chromosome)

        current_path: List[Tuple[str, str]] = []
        last_pos = -1

        for pos in range(len(species_busco)):
            if pos not in gene_chroms:
                if len(current_path) >= 2:
                    observed_chains.add(tuple(current_path))
                current_path = []
                last_pos = -1
                continue

            node = gene_chroms[pos]

            if not current_path:
                current_path = [node]
                last_pos = pos
                continue

            if pos == last_pos + 1 and (gene_chroms[last_pos], node) in sig_edges:
                current_path.append(node)
                last_pos = pos
            else:
                if len(current_path) >= 2:
                    observed_chains.add(tuple(current_path))
                current_path = [node]
                last_pos = pos

        if len(current_path) >= 2:
            observed_chains.add(tuple(current_path))

    # Remove sub-chains that are contiguous sub-paths of longer chains.
    # A gene walk can emit a fragment (e.g. [(A,chr1),(B,chr2)]) that is a
    # prefix/suffix/interior of a longer chain produced by another gene.
    # Keeping both would let genes match the shorter chain even when they
    # diverge from the longer one at positions the short chain doesn't cover.
    filtered: Set[Tuple[Tuple[str, str], ...]] = set()
    for chain in observed_chains:
        is_subchain = False
        for other in observed_chains:
            if len(other) <= len(chain):
                continue
            # Check if chain appears as a contiguous subsequence of other
            for start in range(len(other) - len(chain) + 1):
                if other[start : start + len(chain)] == chain:
                    is_subchain = True
                    break
            if is_subchain:
                break
        if not is_subchain:
            filtered.add(chain)

    chains = [list(chain) for chain in sorted(filtered)]
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

    eligible = _eligible_genes(species_busco)

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
