import argparse
import glob
import os
from collections import defaultdict
from typing import Dict, List, Tuple

from hobrac.rename_chr import fasta_basename

from .coloring import apply_custom_colors, apply_custom_colors_with_algs
from .io import (
    parse_custom_algs,
    parse_custom_colors,
    parse_organism_name,
    read_busco_tsv,
    read_fasta_sizes,
)
from .models import (
    DEFAULT_COLOR,
    BuscoGene,
    ChromosomeAssociation,
)
from .ordering import (
    get_chromosome_order,
    get_chromosome_order_by_gravity,
    get_chromosome_order_by_span,
)
from .rearrangement import calculate_rearrangement_indices, save_rearrangement_indices
from .statistics import detect_algs_transitive


def generate_bed_file(
    busco_data: Dict[str, BuscoGene], species_name: str, output_path: str
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
        busco_data.values(), key=lambda g: (g.chromosome, min(g.start, g.end))
    )

    with open(output_path, "w") as f:
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
    hide_non_significant: bool = False,
) -> None:
    """
    Generate a JCVI .simple file with one line per BUSCO gene pair.

    Each line represents a single ortholog and uses its ``gene_colors`` entry
    as the link color. Because ``gene_colors`` is derived from the single,
    run-wide ``gene_to_chain`` mapping, this guarantees that a gene keeps the
    same color across every .simple file of a single run — no block-level
    majority vote can override its assignment.

    JCVI karyotype expects the 6-column format:
    ``startGene1 endGene1 startGene2 endGene2 score orientation``

    Args:
        busco1: BUSCO data for species 1
        busco2: BUSCO data for species 2
        gene_colors: Dictionary mapping BUSCO ID to color
        sp1: Species 1 name
        sp2: Species 2 name
        output_path: Output .simple file path
        hide_non_significant: If True, skip genes colored lightgrey
    """
    common_ids = set(busco1.keys()) & set(busco2.keys())
    if not common_ids:
        open(output_path, "w").close()
        return

    # Non-significant (grey) genes are written first so that JCVI draws the
    # significant ribbons on top of them. Ties are broken by gene id for
    # deterministic output.
    sorted_ids = sorted(
        common_ids,
        key=lambda gid: (
            gene_colors.get(gid, DEFAULT_COLOR) != DEFAULT_COLOR,
            gid,
        ),
    )

    with open(output_path, "w") as f:
        for busco_id in sorted_ids:
            color = gene_colors.get(busco_id, DEFAULT_COLOR)
            if hide_non_significant and color == DEFAULT_COLOR:
                continue
            gene1 = f"{sp1}_{busco_id}"
            gene2 = f"{sp2}_{busco_id}"
            f.write(f"{color}*{gene1}\t{gene1}\t{gene2}\t{gene2}\t1\t+\n")


def generate_seqids_file(
    species_data: List[Tuple[str, Dict[str, BuscoGene]]],
    output_path: str,
    use_gravity_ordering: bool = False,
    assembly_fasta_sizes: Dict[str, int] = None,
) -> List[List[str]]:
    """
    Generate JCVI seqids file with chromosome order for each species.

    Args:
        species_data: List of (species_name, busco_data) tuples
        output_path: Output seqids file path
        use_gravity_ordering: If True, order subsequent species by gravity
                              with the species above them
        assembly_fasta_sizes: Dictionary mapping assembly sequence names to sizes.
                              Used for ordering the first species (assembly).

    Returns:
        The per-species chromosome order (one list per species, in the same
        order as ``species_data``). This is the exact axis order the ALG
        dotplots reuse so they line up with the karyotype.
    """
    all_orders: List[List[str]] = []
    with open(output_path, "w") as f:
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
            all_orders.append(chromosomes)

            # Store for next iteration
            previous_order = chromosomes
            previous_busco = busco_data

    return all_orders


def save_dotplot_axis_orders(
    species_keys: List[str],
    orders: List[List[str]],
    output_dir: str,
) -> None:
    """
    Write one axis-order file per species for the ALG dotplots.

    Each file is named ``<key>.order`` and holds the species' comma-separated
    chromosome order — identical to its ``seqids`` line. ``species_keys`` are the
    stable identifiers (assembly name / reference accession) captured *before*
    any custom ``--jcvi-names`` renaming, so the dotplot rule can look up a
    species by its accession wildcard regardless of display name.
    """
    orders_dir = os.path.join(output_dir, "dotplot_orders")
    os.makedirs(orders_dir, exist_ok=True)
    for key, order in zip(species_keys, orders):
        with open(os.path.join(orders_dir, f"{key}.order"), "w") as f:
            f.write(",".join(order) + "\n")


def resolve_layout_labels(
    species_keys: List[str],
    assembly_display_name: str,
    custom_names: str,
    reference_dir: str,
) -> List[str]:
    """Resolve one karyotype track label per species, mirroring the dotplot grid.

    Kept in lockstep with ``grid.resolve_titles`` / ``grid.resolve_assembly_title``
    so the karyotype and the dotplot grid show the same names. Priority per
    species (first hit wins):

      1. ``custom_names`` (``--jcvi-names``): a full ordered list
         (1 assembly + N references) names every track; a single name targets the
         assembly only, leaving references to fall through.
      2. assembly -> ``assembly_display_name`` (hobrac's ``-n/--name``, spaces
         preserved); references -> the NCBI ``# Organism name:`` from
         ``reference/<accession>_assembly_report.txt``.
      3. references -> the accession itself.

    ``species_keys`` are the stable, underscore-free keys (underscored assembly
    name, then reference accessions). Labels are display-only: they never feed the
    bed/links gene-id prefixes, so the assembly label can carry spaces while its
    key stays underscored.
    """
    names_list = [n.strip() for n in custom_names.split(",")] if custom_names else []
    full_list = len(names_list) == len(species_keys)

    labels: List[str] = []
    for i, key in enumerate(species_keys):
        if full_list and names_list[i]:
            labels.append(names_list[i])
        elif i == 0:
            # Assembly: a single --jcvi-names entry targets it; else spaced name.
            labels.append(
                names_list[0] if names_list and names_list[0] else assembly_display_name
            )
        else:
            report = os.path.join(reference_dir, f"{key}_assembly_report.txt")
            labels.append(parse_organism_name(report) or key)
    return labels


def generate_layouts_file(
    species_beds: List[Tuple[str, str]],
    links_files: List[str],
    output_path: str,
    labels: List[str] = None,
) -> None:
    """
    Generate JCVI layouts file.

    Args:
        species_beds: List of (species_name, bed_file_path) tuples
        links_files: List of links file paths
        output_path: Output layouts file path
        labels: Display label per track, in ``species_beds`` order. Falls back to
                ``species_name`` when omitted. jcvi's own ``label`` column is left
                blank so long names cannot overlap the chromosome tracks; the
                labels are stashed in a jcvi-ignored ``# track_labels`` comment
                instead, and ``karyotype_legend`` draws them right-aligned (wrapped
                to two lines when needed) as a post-processing step.
    """
    num_species = len(species_beds)
    track_labels = labels if labels else [name for name, _ in species_beds]

    with open(output_path, "w") as f:
        f.write("# y, xstart, xend, rotation, color, label, va, bed\n")
        f.write(f"#{'-' * 60}\n")
        # jcvi ignores '#' lines; the post-processor reads the labels from here.
        f.write("# track_labels\t" + "\t".join(track_labels) + "\n")

        # Calculate y positions (evenly spaced from 0.9 to 0.1)
        y_positions = [
            0.9 - (i * 0.8 / max(num_species - 1, 1)) for i in range(num_species)
        ]

        for i, (species_name, bed_path) in enumerate(species_beds):
            bed_basename = os.path.basename(bed_path)
            y = y_positions[i]
            va = "top" if i == 0 else "bottom"
            # Empty label column: labels are drawn by karyotype_legend instead.
            f.write(
                f"{y:.2f},\t0.1,\t0.96,\t0,\tblack,\t,\t{va},\t{bed_basename}\n"
            )

        f.write("\n# edges\n")
        for i, links_file in enumerate(links_files):
            links_basename = os.path.basename(links_file)
            f.write(f"e, {i}, {i + 1}, {links_basename}\n")


def save_chromosome_associations(
    associations: List[ChromosomeAssociation], output_path: str
) -> None:
    """
    Save chromosome associations to a TSV file with header.

    Args:
        associations: List of chromosome associations
        output_path: Output TSV file path
    """
    with open(output_path, "w") as f:
        f.write(
            "species1\tspecies2\tchr1\tchr2\tp_value\t"
            "corrected_p_value\tgene_count\tsignificant\n"
        )
        for a in associations:
            significant = "accepted" if a.significant else "rejected"
            f.write(
                f"{a.species1}\t{a.species2}\t{a.chr1}\t{a.chr2}\t"
                f"{a.p_value:.2e}\t{a.corrected_p_value:.2e}\t"
                f"{a.gene_count}\t{significant}\n"
            )


def save_chains(
    chains: List[List[Tuple[str, str]]],
    chain_colors: Dict[int, str],
    gene_to_chain: Dict[str, int],
    species_names: List[str],
    output_path: str,
) -> None:
    """
    Save chromosome chains to a TSV file in wide format.

    One row per chain. Columns: chain_id, color (palette hex), n_genes,
    then one column per species containing the chromosome covered by that
    chain, or '-' when the chain does not cover that species. Species
    columns appear in *species_names* order.

    The file is always written, even when *chains* is empty — in that case
    only the header row is emitted.
    """
    gene_counts: Dict[int, int] = defaultdict(int)
    for chain_id in gene_to_chain.values():
        if chain_id >= 0:
            gene_counts[chain_id] += 1

    header = ["chain_id", "color", "n_genes", *species_names]
    with open(output_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for chain_id, chain in enumerate(chains):
            chrom_by_species = {sp: chrom for sp, chrom in chain}
            row = [
                str(chain_id),
                chain_colors[chain_id],
                str(gene_counts[chain_id]),
                *(chrom_by_species.get(sp, "-") for sp in species_names),
            ]
            f.write("\t".join(row) + "\n")


def save_gene_chains(
    species_busco: List[Tuple[str, Dict[str, BuscoGene]]],
    gene_to_chain: Dict[str, int],
    gene_colors: Dict[str, str],
    custom_algs: Dict[str, str],
    output_path: str,
) -> None:
    """
    Save a per-gene chain table to a TSV file in wide format.

    One row per BUSCO gene present in at least one species. Columns:
    ``chain_id``, ``custom_alg_id``, ``gene``, ``color`` (the gene's run-wide
    identity color), one column per species (in *species_busco* order) holding
    the gene's own chromosome in that species or ``ABSENT`` when the gene is
    missing there, then one ``<species>_pos`` column per species holding the
    gene's coordinates as ``start:end`` (raw BUSCO ints), or ``ABSENT`` when
    missing.

    ``chain_id`` is hobrac's own chain index, or ``-`` when the gene is on no
    chain. ``custom_alg_id`` is the reference ALG label supplied by the user in
    the custom color file (its third column), or ``-`` when there is no color
    file or the gene is not listed with an ALG. The two sit side by side so the
    hobrac-detected chain can be compared against the reference ALG.

    ``gene_colors`` is the run-wide identity color resolved the same way the
    karyotype renders it (custom override, else chain palette hex, else
    ``lightgrey``); genes missing from it fall back to ``lightgrey``. Because it
    is run-wide it cannot encode the plot's per-pair ``chains_covering_pair``
    gating: a gene whose chain spans only some species pairs keeps its color
    here but is drawn grey on the pairs its chain does not cover. See
    ``resolve_gene_identity_colors``.

    Rows are ordered with chain genes first, grouped by ``chain_id`` then gene
    id, and no-chain (grey) genes last by gene id. The file is always written,
    even when no genes are present — in that case only the header is emitted.
    """
    custom_algs = custom_algs or {}
    species_names = [name for name, _ in species_busco]

    all_ids = set()
    for _, busco_data in species_busco:
        all_ids.update(busco_data.keys())

    def sort_key(gid: str) -> Tuple[bool, int, str]:
        chain_id = gene_to_chain.get(gid, -1)
        return (chain_id < 0, chain_id if chain_id >= 0 else 0, gid)

    pos_names = [f"{name}_pos" for name in species_names]
    header = [
        "chain_id",
        "custom_alg_id",
        "gene",
        "color",
        *species_names,
        *pos_names,
    ]
    with open(output_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for gid in sorted(all_ids, key=sort_key):
            chain_id = gene_to_chain.get(gid, -1)
            chroms = [
                busco_data[gid].chromosome if gid in busco_data else "ABSENT"
                for _, busco_data in species_busco
            ]
            positions = [
                (
                    f"{busco_data[gid].start}:{busco_data[gid].end}"
                    if gid in busco_data
                    else "ABSENT"
                )
                for _, busco_data in species_busco
            ]
            custom_alg = custom_algs.get(gid, "-")
            color = gene_colors.get(gid, DEFAULT_COLOR)
            chain_label = str(chain_id) if chain_id >= 0 else "-"
            row = [chain_label, custom_alg, gid, color, *chroms, *positions]
            f.write("\t".join(row) + "\n")


def resolve_gene_identity_colors(
    all_gene_ids: set,
    gene_to_chain: Dict[str, int],
    chain_colors: Dict[int, str],
    custom_colors: Dict[str, str],
    skip_alg: bool,
) -> Dict[str, str]:
    """
    Resolve one run-wide identity color per gene, mirroring the plot.

    This is the run-wide counterpart of ``apply_custom_colors_with_algs``: it
    resolves a single color per gene using the same fallback rules the
    karyotype applies, minus the per-pair ``chains_covering_pair`` gating that a
    one-row-per-gene table cannot represent.

    Rules:
      * ``custom_colors`` provided + ``skip_alg`` (no chains): custom color, or
        ``lightgrey`` when the gene is not listed.
      * gene on no chain (``chain_id < 0``): ``lightgrey``.
      * gene on a chain, ``custom_colors`` provided: custom color, or
        ``lightgrey`` when the gene is not listed — matching the plot, which
        never falls back to the chain palette once a custom file is given.
      * gene on a chain, no ``custom_colors``: the chain palette hex.

    Note the residual run-wide vs per-pair difference: a gene whose chain covers
    only some species pairs keeps its custom color here but is drawn grey on the
    pairs its chain does not span.

    Args:
        all_gene_ids: Every BUSCO id present in at least one species.
        gene_to_chain: Dict mapping BUSCO gene id to chain id (-1 = no chain).
        chain_colors: Dict mapping chain id to hex color.
        custom_colors: Dictionary mapping BUSCO id to hex color (may be empty).
        skip_alg: Whether ALG/chain detection was skipped.

    Returns:
        Dictionary mapping BUSCO id to its run-wide identity color.
    """
    gene_identity_colors: Dict[str, str] = {}
    for gid in all_gene_ids:
        if custom_colors and skip_alg:
            gene_identity_colors[gid] = custom_colors.get(gid, DEFAULT_COLOR)
        else:
            chain_id = gene_to_chain.get(gid, -1)
            if chain_id < 0:
                gene_identity_colors[gid] = DEFAULT_COLOR
            elif custom_colors:
                gene_identity_colors[gid] = custom_colors.get(gid, DEFAULT_COLOR)
            else:
                gene_identity_colors[gid] = chain_colors[chain_id]
    return gene_identity_colors


def run(
    assembly_busco: str,
    assembly_fasta: str,
    busco_refs: List[str],
    accession_order: str,
    manual_refs: str,
    output_dir: str,
    assembly_name: str = "assembly",
    assembly_display_name: str = "",
    reference_dir: str = "reference",
    use_gravity_ordering: bool = True,
    min_busco_genes: int = 0,
    custom_color_file: str = "",
    custom_names: str = "",
    hide_non_significant: bool = False,
    skip_alg: bool = False,
    alpha: float = 0.01,
    min_chain_genes: int = 5,
    permissive_alg: bool = False,
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
        assembly_name: Stable, underscore-free key for the assembly. Used as the
                       gene-id prefix in bed/links files and as the dotplot-order
                       key, so it must not contain whitespace.
        assembly_display_name: Human-facing assembly label for the karyotype
                       (hobrac's ``-n/--name``, spaces preserved). Falls back to
                       ``assembly_name`` when empty. Display-only.
        reference_dir: Directory holding ``<accession>_assembly_report.txt`` NCBI
                       reports, used to label references by organism name.
        use_gravity_ordering: If True, order chromosomes by gravity for
                              diagonal alignment patterns
        min_busco_genes: Minimum complete BUSCO genes required per sequence
        custom_color_file: Path to custom color file for gene coloring.
        custom_names: Comma-separated custom names for JCVI tracks. If one name
                      is provided, it applies only to the assembly. If multiple
                      names are provided, the count must match species count.
        hide_non_significant: If True, hide links between chromosome pairs
                              without significant associations.
        skip_alg: If True, skip ALG statistical testing when using custom
                  colors. Only effective with custom_color_file.
        min_chain_genes: Minimum BUSCO genes a chain must be supported by.
                         Chains with fewer supporting genes are dropped before
                         sub-chain pruning.

    Returns:
        Dictionary with paths to generated files
    """
    os.makedirs(output_dir, exist_ok=True)

    assembly_fasta_sizes = read_fasta_sizes(assembly_fasta)

    custom_colors = {}
    custom_algs = {}
    if custom_color_file:
        custom_colors = parse_custom_colors(custom_color_file)
        custom_algs = parse_custom_algs(custom_color_file)

    # Build mapping from accession to BUSCO path
    ref_busco_paths = {}
    for ref_dir in busco_refs:
        # Extract accession from directory name (busco_reference_ACCESSION)
        dirname = os.path.basename(ref_dir.rstrip("/"))
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
                name = fasta_basename(ref_path)
                if name in ref_busco_paths:
                    species_order.append((name, ref_busco_paths[name]))
    else:
        with open(accession_order) as f:
            for line in f:
                accession = line.strip()
                if accession and accession in ref_busco_paths:
                    species_order.append((accession, ref_busco_paths[accession]))

    # Stable, underscore-free per-species keys (assembly name / reference
    # accession). These key every generated file and gene-id prefix; custom
    # --jcvi-names are display-only and applied later as karyotype/dotplot
    # labels, so names with spaces never leak into jcvi's whitespace-delimited
    # bed/links parsers (which would crash with "too many values to unpack").
    species_keys = [name for name, _ in species_order]

    species_busco = []
    for name, path in species_order:
        busco_data = read_busco_tsv(path, min_busco_genes)
        species_busco.append((name, busco_data))

    rearrangement_summary = ""
    rearrangement_by_alg = ""
    if custom_algs:
        rearrangement_summary = os.path.join(output_dir, "rearrangement_index.tsv")
        rearrangement_by_alg = os.path.join(
            output_dir, "rearrangement_index_by_alg.tsv"
        )
        rearrangement_rows = {
            name: calculate_rearrangement_indices(name, busco_data, custom_algs)
            for name, busco_data in species_busco
        }
        save_rearrangement_indices(
            rearrangement_rows,
            rearrangement_summary,
            rearrangement_by_alg,
        )

    bed_files = []
    for name, busco_data in species_busco:
        bed_path = os.path.join(output_dir, f"{name}.bed")
        generate_bed_file(busco_data, name, bed_path)
        bed_files.append((name, bed_path))

    # Detect ALGs and generate links for consecutive species pairs
    associations_output = os.path.join(output_dir, "chromosome_associations.tsv")

    # Phase 1: Transitive ALG detection across all species
    # Skip only when custom colors are provided AND skip_alg is True
    all_chromosome_associations: List[ChromosomeAssociation] = []
    gene_to_chain: Dict[str, int] = {}
    chain_colors: Dict[int, str] = {}
    chains: List[List[Tuple[str, str]]] = []
    if not (custom_colors and skip_alg):
        _, chains, gene_to_chain, chain_colors, all_chromosome_associations = (
            detect_algs_transitive(
                species_busco,
                alpha=alpha,
                min_chain_genes=min_chain_genes,
                permissive_alg=permissive_alg,
            )
        )
    save_chromosome_associations(all_chromosome_associations, associations_output)

    algs_output = os.path.join(output_dir, "algs.tsv")
    save_chains(
        chains,
        chain_colors,
        gene_to_chain,
        [name for name, _ in species_busco],
        algs_output,
    )

    # Resolve one run-wide identity color per gene, matching how the karyotype
    # renders it (see resolve_gene_identity_colors).
    all_gene_ids = set()
    for _, busco_data in species_busco:
        all_gene_ids.update(busco_data.keys())
    gene_identity_colors = resolve_gene_identity_colors(
        all_gene_ids, gene_to_chain, chain_colors, custom_colors, skip_alg
    )

    gene_chains_output = os.path.join(output_dir, "gene_chains.tsv")
    save_gene_chains(
        species_busco,
        gene_to_chain,
        gene_identity_colors,
        custom_algs,
        gene_chains_output,
    )

    # Phase 2: Generate outputs for each pair
    links_files = []
    for i in range(len(species_busco) - 1):
        sp1_name, sp1_busco = species_busco[i]
        sp2_name, sp2_busco = species_busco[i + 1]

        if custom_colors and skip_alg:
            gene_colors = apply_custom_colors(sp1_busco, sp2_busco, custom_colors)
        else:
            gene_colors = apply_custom_colors_with_algs(
                sp1_busco,
                sp2_busco,
                gene_to_chain,
                chain_colors,
                custom_colors,
                chains=chains,
                sp1_name=sp1_name,
                sp2_name=sp2_name,
            )

        links_path = os.path.join(output_dir, f"links.{sp1_name}.{sp2_name}.simple")
        generate_links_file(
            sp1_busco,
            sp2_busco,
            gene_colors,
            sp1_name,
            sp2_name,
            links_path,
            hide_non_significant=hide_non_significant,
        )
        links_files.append(links_path)

    seqids_path = os.path.join(output_dir, "seqids")
    axis_orders = generate_seqids_file(
        species_busco, seqids_path, use_gravity_ordering, assembly_fasta_sizes
    )

    # Per-species axis-order files keyed by stable accession, reused by the ALG
    # dotplots so their axes match the karyotype exactly.
    save_dotplot_axis_orders(species_keys, axis_orders, output_dir)

    # Resolve display labels for the karyotype tracks the same way the dotplot
    # grid does (see resolve_layout_labels), so the two figures agree. These feed
    # only the layouts `label` column; the stable species_keys still key the
    # bed/links files, so the assembly label may carry spaces while its key stays
    # underscored.
    layout_labels = resolve_layout_labels(
        species_keys,
        assembly_display_name or assembly_name,
        custom_names,
        reference_dir,
    )

    layouts_path = os.path.join(output_dir, "layouts")
    generate_layouts_file(bed_files, links_files, layouts_path, labels=layout_labels)

    return {
        "seqids": seqids_path,
        "layouts": layouts_path,
        "chromosome_associations": associations_output,
        "algs": algs_output,
        "gene_chains": gene_chains_output,
        "rearrangement_index": rearrangement_summary,
        "rearrangement_index_by_alg": rearrangement_by_alg,
        "bed_files": [p for _, p in bed_files],
        "links_files": links_files,
    }


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Generate JCVI karyotype files from BUSCO data with ALG detection"
    )
    parser.add_argument(
        "--busco_assembly", required=True, help="Path to assembly BUSCO full_table.tsv"
    )
    parser.add_argument(
        "--assembly-fasta", required=True, help="Path to the assembly fasta file"
    )
    parser.add_argument(
        "--busco_references",
        nargs="+",
        required=True,
        help="Paths to reference BUSCO directories",
    )
    parser.add_argument(
        "--accession_order",
        required=True,
        help="Path to file with accessions in order (MASH distance)",
    )
    parser.add_argument(
        "--manual_refs", default="", help="Semicolon-separated manual reference paths"
    )
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument(
        "--assembly_name",
        default="assembly",
        help="Stable, underscore-free assembly key (gene-id prefix / dotplot-order "
        "key). Must not contain whitespace.",
    )
    parser.add_argument(
        "--assembly-display-name",
        default="",
        help="Human-facing assembly label for the karyotype (spaces allowed). "
        "Falls back to --assembly_name when omitted.",
    )
    parser.add_argument(
        "--reference-dir",
        default="reference",
        help="Directory holding <accession>_assembly_report.txt NCBI reports, used "
        "to label references by organism name (default: reference).",
    )
    parser.add_argument(
        "--min-busco-genes",
        type=int,
        default=0,
        help="Minimum complete BUSCO genes required per sequence (default: 0)",
    )
    parser.add_argument(
        "--jcvi-custom-colors",
        default="",
        help="Path to custom color file (tab-separated: BUSCO_ID, COLOR, ALG_NAME; "
        "COLOR may be R,G,B, #rrggbb, or rrggbb). "
        "By default, ALG statistical testing still runs to determine significance; "
        "use --skip-alg to disable it. Genes not in the file will be shown in grey. "
        "The ALG_NAME column is also used for rearrangement "
        "index calculation.",
    )
    parser.add_argument(
        "--jcvi-names",
        default="",
        help="Comma-separated custom names for JCVI tracks. If one name is provided, "
        "it applies only to the assembly. If multiple names are provided, the count "
        "must equal 1 (assembly) + number of references, in order.",
    )
    parser.add_argument(
        "--hide-non-significant",
        action="store_true",
        help="Hide links between chromosome pairs without significant associations. "
        "This produces a cleaner plot showing only ALG-related synteny.",
    )
    parser.add_argument(
        "--skip-alg",
        action="store_true",
        help="Skip ALG statistical testing when using custom colors. All genes in the "
        "color file get their custom color; unlisted genes are grey. Only effective "
        "with --jcvi-custom-colors.",
    )
    parser.add_argument(
        "--alg-pvalue",
        type=float,
        default=0.01,
        help="Significance threshold (alpha) for Fisher's exact test in ALG "
        "detection (default: 0.01).",
    )
    parser.add_argument(
        "--jcvi-min-chain-genes",
        type=int,
        default=5,
        help="Minimum BUSCO genes a chromosome chain must be supported by. "
        "Chains with fewer genes are dropped before sub-chain pruning, so "
        "sub-chains hidden inside a dropped long chain can re-emerge "
        "(default: 5).",
    )
    parser.add_argument(
        "--jcvi-permissive-alg",
        action="store_true",
        default=False,
        help="Use a permissive threshold for chain validation. By default each "
        "node in a chain of length n must have significant associations "
        "with at least n/2 other nodes; this flag relaxes the requirement "
        "to just 1 significant link per node.",
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
        assembly_display_name=args.assembly_display_name,
        reference_dir=args.reference_dir,
        min_busco_genes=args.min_busco_genes,
        custom_color_file=args.jcvi_custom_colors,
        custom_names=args.jcvi_names,
        hide_non_significant=args.hide_non_significant,
        skip_alg=args.skip_alg,
        alpha=args.alg_pvalue,
        min_chain_genes=args.jcvi_min_chain_genes,
        permissive_alg=args.jcvi_permissive_alg,
    )


if __name__ == "__main__":
    main()
