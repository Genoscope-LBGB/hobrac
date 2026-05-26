from collections import defaultdict
from typing import Dict, List, Set, Tuple

from .models import BuscoGene, ChromosomeAssociation, PairwiseAssociation


def _eligible_genes(
    species_busco: List[Tuple[str, Dict[str, BuscoGene]]],
) -> Set[str]:
    """Return gene IDs present in at least two consecutive species."""
    eligible: Set[str] = set()
    for i in range(len(species_busco) - 1):
        eligible |= species_busco[i][1].keys() & species_busco[i + 1][1].keys()
    return eligible


def _walk_chain_counts(
    pairwise_associations: List[PairwiseAssociation],
    species_busco: List[Tuple[str, Dict[str, BuscoGene]]],
) -> Dict[Tuple[Tuple[str, str], ...], int]:
    """
    Walk each eligible gene through species and count maximal contiguous
    paths of significant edges. Return a mapping from chain tuple to the
    number of distinct genes whose maximal walk produced that exact chain.

    A gene contributes to a chain only when the chain equals (not merely
    contains) one of its maximal walks. Genes whose walks are strictly
    shorter do not contribute to longer chains' counts.
    """
    if not pairwise_associations or not species_busco:
        return {}

    sig_edges: Set[Tuple[Tuple[str, str], Tuple[str, str]]] = set()
    for assoc in pairwise_associations:
        sig_edges.add(((assoc.species1, assoc.chr1), (assoc.species2, assoc.chr2)))

    eligible = _eligible_genes(species_busco)
    chain_counts: Dict[Tuple[Tuple[str, str], ...], int] = defaultdict(int)

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
                    chain_counts[tuple(current_path)] += 1
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
                    chain_counts[tuple(current_path)] += 1
                current_path = [node]
                last_pos = pos

        if len(current_path) >= 2:
            chain_counts[tuple(current_path)] += 1

    return chain_counts


def enumerate_chains(
    pairwise_associations: List[PairwiseAssociation],
    species_busco: List[Tuple[str, Dict[str, BuscoGene]]],
    min_chain_genes: int = 0,
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
        min_chain_genes: Drop chains whose maximal-walk support is strictly
            below this threshold before sub-chain pruning. Support is the
            number of distinct genes whose maximal walk equals the chain
            exactly — distinct from the partial-match count used by
            :func:`build_gene_chain_mapping`. Default 0 (keep every chain).

    Returns:
        List of chains sorted deterministically, where each chain is a list
        of (species, chromosome) tuples in species order.
    """
    chain_counts = _walk_chain_counts(pairwise_associations, species_busco)
    surviving = {
        chain for chain, count in chain_counts.items() if count >= min_chain_genes
    }
    surviving = _prune_subchains(surviving)
    return [list(chain) for chain in sorted(surviving)]


def _prune_subchains(
    chains: Set[Tuple[Tuple[str, str], ...]],
) -> Set[Tuple[Tuple[str, str], ...]]:
    """
    Remove chains that are contiguous sub-paths of longer chains.

    A gene walk can emit a fragment (e.g. [(A,chr1),(B,chr2)]) that is a
    prefix/suffix/interior of a longer chain produced by another gene.
    Keeping both would let genes match the shorter chain even when they
    diverge from the longer one at positions the short chain doesn't cover.
    """
    filtered: Set[Tuple[Tuple[str, str], ...]] = set()
    for chain in chains:
        is_subchain = False
        for other in chains:
            if len(other) <= len(chain):
                continue
            for start in range(len(other) - len(chain) + 1):
                if other[start : start + len(chain)] == chain:
                    is_subchain = True
                    break
            if is_subchain:
                break
        if not is_subchain:
            filtered.add(chain)
    return filtered


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


def _count_matching_genes(
    chain: List[Tuple[str, str]],
    species_busco: List[Tuple[str, Dict[str, BuscoGene]]],
    species_index: Dict[str, int],
    eligible: Set[str],
) -> int:
    chain_path = {species_index[sp]: chrom for sp, chrom in chain}
    count = 0
    for gene_id in eligible:
        if _gene_matches_chain(
            {
                pos: busco_data[gene_id].chromosome
                for pos, (_, busco_data) in enumerate(species_busco)
                if gene_id in busco_data
            },
            chain_path,
        ):
            count += 1
    return count


def _count_sig_links(
    node: Tuple[str, str],
    nodes: List[Tuple[str, str]],
    sig_pairs: Set[Tuple[Tuple[str, str], Tuple[str, str]]],
) -> int:
    return sum(1 for other in nodes if other is not node and (node, other) in sig_pairs)


def validate_chains(
    chains: List[List[Tuple[str, str]]],
    all_chromosome_associations: List[ChromosomeAssociation],
    permissive: bool = False,
    min_chain_genes: int = 0,
    species_busco: List[Tuple[str, Dict[str, BuscoGene]]] = None,
) -> List[List[Tuple[str, str]]]:
    """
    Validate chains by iteratively pruning nodes that lack enough significant
    links to other surviving nodes. After pruning, remaining contiguous
    segments of length >= 2 are kept and deduplicated.

    Args:
        chains: Chains from enumerate_chains
        all_chromosome_associations: All tested ChromosomeAssociation objects
            (significant and non-significant) across all species pairs
        permissive: If True, require only 1 significant link per node.
            If False (default), require int(n/2) where n is the chain length.
        min_chain_genes: Drop sub-chains whose gene support is below this
            threshold (re-checked after splitting).
        species_busco: Needed when min_chain_genes > 0 to recount gene
            support on sub-chains.

    Returns:
        Validated (and possibly split) chains, sorted deterministically.
    """
    sig_pairs: Set[Tuple[Tuple[str, str], Tuple[str, str]]] = set()
    for assoc in all_chromosome_associations:
        if assoc.significant:
            a = (assoc.species1, assoc.chr1)
            b = (assoc.species2, assoc.chr2)
            sig_pairs.add((a, b))
            sig_pairs.add((b, a))

    validated: List[List[Tuple[str, str]]] = []
    for chain in chains:
        surviving = list(chain)
        while True:
            n = len(surviving)
            if n < 2:
                break
            threshold = 1 if permissive else n // 2
            still_ok = [
                node
                for node in surviving
                if _count_sig_links(node, surviving, sig_pairs) >= threshold
            ]
            if len(still_ok) == len(surviving):
                break
            surviving = still_ok

        if len(surviving) < 2:
            continue

        original_indices = [chain.index(node) for node in surviving]
        segment: List[Tuple[str, str]] = []
        prev = -2
        for k, idx in enumerate(original_indices):
            if idx == prev + 1:
                segment.append(surviving[k])
            else:
                if len(segment) >= 2:
                    validated.append(segment)
                segment = [surviving[k]]
            prev = idx
        if len(segment) >= 2:
            validated.append(segment)

    seen: Set[Tuple[Tuple[str, str], ...]] = set()
    deduped: List[List[Tuple[str, str]]] = []
    for seg in validated:
        key = tuple(seg)
        if key not in seen:
            seen.add(key)
            deduped.append(seg)
    validated = deduped

    if min_chain_genes > 0 and species_busco:
        species_index = {name: i for i, (name, _) in enumerate(species_busco)}
        eligible = _eligible_genes(species_busco)
        validated = [
            seg
            for seg in validated
            if _count_matching_genes(seg, species_busco, species_index, eligible)
            >= min_chain_genes
        ]

    validated = [list(c) for c in _prune_subchains({tuple(c) for c in validated})]
    return sorted(validated)
