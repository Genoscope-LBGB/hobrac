from collections import defaultdict
from typing import Dict, List, Set, Tuple

from .models import BuscoGene, ChromosomeAssociation, PairwiseAssociation

Node = Tuple[str, str]
Chain = Tuple[Node, ...]
ChainEvidence = Dict[Chain, Set[str]]

# Refuse incomplete, name-dependent results if exact enumeration is too costly.
_MAX_DFS_VISITS = 500_000


def _eligible_genes(
    species_busco: List[Tuple[str, Dict[str, BuscoGene]]],
) -> Set[str]:
    """Return gene IDs present in at least two species."""
    counts: Dict[str, int] = defaultdict(int)
    for _, busco_data in species_busco:
        for gene_id in busco_data:
            counts[gene_id] += 1
    return {g for g, c in counts.items() if c >= 2}


def _longest_paths_per_component(
    nodes: List[Node], adj: Dict[Node, Set[Node]]
) -> List[Chain]:
    """Return every distinct longest-path node set in each induced component.

    Node order is only a canonical identifier, not a selected traversal. Search
    states merge traversals with the same visited nodes and endpoint. Once a
    path covers a whole component, its node set is the sole possible result.
    """
    remaining = set(nodes)
    result: List[Chain] = []
    while remaining:
        root = min(remaining)
        remaining.remove(root)
        stack = [root]
        component = []
        while stack:
            node = stack.pop()
            component.append(node)
            neighbors = adj.get(node, set()) & remaining
            remaining.difference_update(neighbors)
            stack.extend(neighbors)

        if len(component) < 2:
            continue
        component.sort()
        indices = {node: i for i, node in enumerate(component)}
        neighbor_masks = [
            sum(
                1 << indices[neighbor]
                for neighbor in adj.get(node, set())
                if neighbor in indices
            )
            for node in component
        ]
        full_mask = (1 << len(component)) - 1
        search = [(i, 1 << i) for i in range(len(component))]
        seen = set(search)
        best_length = 0
        best_masks: Set[int] = set()
        visits = 0
        while search:
            endpoint, mask = search.pop()
            visits += 1
            if visits > _MAX_DFS_VISITS:
                raise RuntimeError(
                    f"Exact ALG path search exceeded {_MAX_DFS_VISITS} states "
                    f"for a component of {len(component)} chromosomes; "
                    "refusing to return incomplete chain support."
                )
            length = mask.bit_count()
            if length > best_length:
                best_length = length
                best_masks = {mask}
            elif length == best_length:
                best_masks.add(mask)
            if mask == full_mask:
                break
            neighbors = neighbor_masks[endpoint] & ~mask
            while neighbors:
                bit = neighbors & -neighbors
                neighbors ^= bit
                state = (bit.bit_length() - 1, mask | bit)
                if state not in seen:
                    seen.add(state)
                    search.append(state)
        result.extend(
            tuple(node for i, node in enumerate(component) if mask & (1 << i))
            for mask in best_masks
        )
    return sorted(result)


def enumerate_chains(
    pairwise_associations: List[PairwiseAssociation],
    species_busco: List[Tuple[str, Dict[str, BuscoGene]]],
    min_chain_genes: int = 0,
) -> ChainEvidence:
    """Discover chains and their distinct observed-path supporting genes.

    Each gene supports every longest-path node set in its significant-link
    graph, once per set regardless of traversal count. Missing species never
    supply support. Thresholding precedes removal of strict node-set subsets,
    so a shorter chain can survive when a weak longer chain is filtered out.
    """
    adj: Dict[Node, Set[Node]] = defaultdict(set)
    for assoc in pairwise_associations:
        a = (assoc.species1, assoc.chr1)
        b = (assoc.species2, assoc.chr2)
        adj[a].add(b)
        adj[b].add(a)

    # Genes at the same chromosome locations have the same induced graph.
    locations: ChainEvidence = defaultdict(set)
    for gene_id in _eligible_genes(species_busco):
        nodes = tuple(
            sorted(
                (species, data[gene_id].chromosome)
                for species, data in species_busco
                if gene_id in data
            )
        )
        locations[nodes].add(gene_id)

    evidence: ChainEvidence = defaultdict(set)
    for nodes, genes in locations.items():
        for chain in _longest_paths_per_component(list(nodes), adj):
            evidence[chain].update(genes)
    return _prune_subchains(
        {
            chain: genes
            for chain, genes in evidence.items()
            if len(genes) >= min_chain_genes
        }
    )


def _prune_subchains(evidence: ChainEvidence) -> ChainEvidence:
    """Remove strict node-set subsets, without transferring their support."""
    node_sets = {chain: frozenset(chain) for chain in evidence}
    return {
        chain: evidence[chain]
        for chain in sorted(evidence)
        if not any(
            len(other) > len(chain) and node_sets[chain] < node_sets[other]
            for other in evidence
        )
    }


def build_gene_chain_mapping(
    species_busco: List[Tuple[str, Dict[str, BuscoGene]]],
    chains: List[List[Node]],
    chain_support: Dict[int, int],
) -> Tuple[Dict[str, List[int]], Dict[str, int]]:
    """Return compatible chains and one display choice per eligible gene.

    Missing species are tolerated for compatibility, not for observed support.
    Rank compatible chains by frozen support, then the number of consecutive
    layout links present for this gene, then chain ID. Choosing a display chain
    never changes support or removes alternative memberships.
    """
    species_index = {name: i for i, (name, _) in enumerate(species_busco)}
    chain_paths = [
        {species_index[sp]: chrom for sp, chrom in chain} for chain in chains
    ]
    chain_links = [
        [(pos, pos + 1) for pos in path if pos + 1 in path] for path in chain_paths
    ]
    gene_to_chains: Dict[str, List[int]] = {}
    gene_to_display_chain: Dict[str, int] = {}
    for gene_id in sorted(_eligible_genes(species_busco)):
        gene_chroms = {
            pos: data[gene_id].chromosome
            for pos, (_, data) in enumerate(species_busco)
            if gene_id in data
        }
        matches = [
            cid
            for cid, path in enumerate(chain_paths)
            if _gene_matches_chain(gene_chroms, path)
        ]
        gene_to_chains[gene_id] = matches
        gene_to_display_chain[gene_id] = min(
            matches,
            key=lambda cid: (
                -chain_support[cid],
                -sum(
                    a in gene_chroms and b in gene_chroms for a, b in chain_links[cid]
                ),
                cid,
            ),
            default=-1,
        )
    return gene_to_chains, gene_to_display_chain


def _gene_matches_chain(
    gene_chroms: Dict[int, str], chain_path: Dict[int, str]
) -> bool:
    """Return True if the gene is compatible at every shared species."""
    compared = False
    for pos, chrom in chain_path.items():
        if pos in gene_chroms:
            if gene_chroms[pos] != chrom:
                return False
            compared = True
    return compared


def validate_chains(
    evidence: ChainEvidence,
    all_chromosome_associations: List[ChromosomeAssociation],
    permissive: bool = False,
    min_chain_genes: int = 0,
) -> ChainEvidence:
    """Prune weak nodes and retain observed longest paths in surviving graphs.

    Strict mode requires floor(n/2) significant links per node; permissive mode
    requires one. Recompute until stable, then split using actual significant
    edges rather than an arbitrary traversal order. Each resulting chain
    inherits its parent candidates' distinct observed supporters; unioning
    prevents a shared gene from being counted twice when candidates converge.
    """
    adj: Dict[Node, Set[Node]] = defaultdict(set)
    for assoc in all_chromosome_associations:
        if assoc.significant:
            a = (assoc.species1, assoc.chr1)
            b = (assoc.species2, assoc.chr2)
            adj[a].add(b)
            adj[b].add(a)

    validated: ChainEvidence = defaultdict(set)
    for chain, genes in evidence.items():
        surviving = set(chain)
        while len(surviving) >= 2:
            threshold = 1 if permissive else len(surviving) // 2
            still_ok = {
                node for node in surviving if len(adj[node] & surviving) >= threshold
            }
            if still_ok == surviving:
                break
            surviving = still_ok
        if len(surviving) < 2:
            continue
        for segment in _longest_paths_per_component(list(surviving), adj):
            validated[segment].update(genes)

    return _prune_subchains(
        {
            chain: genes
            for chain, genes in validated.items()
            if len(genes) >= min_chain_genes
        }
    )
