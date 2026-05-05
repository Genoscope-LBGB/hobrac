from collections import defaultdict
from typing import Dict, Iterable, List, Tuple

from .models import BuscoGene, RearrangementAlgIndex


def calculate_rearrangement_indices(
    species: str,
    busco_data: Dict[str, BuscoGene],
    gene_to_alg: Dict[str, str],
) -> List[RearrangementAlgIndex]:
    """
    Calculate rearrangement index components for one species
    (https://academic.oup.com/mbe/article/41/9/msae172/7733614#483558069).

    Only BUSCO genes that are both present in the species and annotated with an
    ALG are used. ALGs absent from a species are not counted in the average,
    because absence reflects marker recovery rather than genome rearrangement.
    """
    chrom_alg_counts: Dict[Tuple[str, str], int] = defaultdict(int)
    alg_counts: Dict[str, int] = defaultdict(int)
    chrom_counts: Dict[str, int] = defaultdict(int)

    for busco_id, gene in busco_data.items():
        alg = gene_to_alg.get(busco_id)
        if not alg:
            continue
        chrom_alg_counts[(gene.chromosome, alg)] += 1
        alg_counts[alg] += 1
        chrom_counts[gene.chromosome] += 1

    rows: List[RearrangementAlgIndex] = []
    for alg in sorted(alg_counts):
        candidates = [
            (chrom, count)
            for (chrom, candidate_alg), count in chrom_alg_counts.items()
            if candidate_alg == alg
        ]
        best_chromosome, best_count = sorted(
            candidates, key=lambda item: (-item[1], item[0])
        )[0]

        splitting_parameter = best_count / alg_counts[alg]
        combining_parameter = best_count / chrom_counts[best_chromosome]
        rearrangement_index = 1 - (splitting_parameter * combining_parameter)

        rows.append(
            RearrangementAlgIndex(
                species=species,
                alg=alg,
                best_chromosome=best_chromosome,
                alg_gene_count=alg_counts[alg],
                chromosome_gene_count=chrom_counts[best_chromosome],
                alg_genes_on_best_chromosome=best_count,
                splitting_parameter=splitting_parameter,
                combining_parameter=combining_parameter,
                rearrangement_index=rearrangement_index,
            )
        )

    return rows


def summarize_rearrangement_indices(
    rows: Iterable[RearrangementAlgIndex],
) -> Tuple[float, float, float, int, int]:
    """Return Ri, splitting index, combining index, n ALGs, and n genes."""
    row_list = list(rows)
    if not row_list:
        return 0.0, 0.0, 0.0, 0, 0

    n_algs = len(row_list)
    n_genes = sum(row.alg_gene_count for row in row_list)
    rearrangement_index = sum(row.rearrangement_index for row in row_list) / n_algs
    splitting_index = sum(1 - row.splitting_parameter for row in row_list) / n_algs
    combining_index = sum(1 - row.combining_parameter for row in row_list) / n_algs
    return rearrangement_index, splitting_index, combining_index, n_algs, n_genes


def save_rearrangement_indices(
    species_rows: Dict[str, List[RearrangementAlgIndex]],
    summary_path: str,
    by_alg_path: str,
) -> None:
    """Save rearrangement index summary and per-ALG detail TSV files."""
    with open(summary_path, "w") as f:
        f.write(
            "species\trearrangement_index\tsplitting_index\tcombining_index\t"
            "n_algs\tn_genes\n"
        )
        for species in sorted(species_rows):
            ri, si, ci, n_algs, n_genes = summarize_rearrangement_indices(
                species_rows[species]
            )
            f.write(
                f"{species}\t{ri:.6f}\t{si:.6f}\t{ci:.6f}\t" f"{n_algs}\t{n_genes}\n"
            )

    with open(by_alg_path, "w") as f:
        f.write(
            "species\talg\tbest_chromosome\talg_gene_count\t"
            "chromosome_gene_count\talg_genes_on_best_chromosome\t"
            "splitting_parameter\tcombining_parameter\t"
            "splitting_index\tcombining_index\trearrangement_index\n"
        )
        for species in sorted(species_rows):
            for row in species_rows[species]:
                f.write(
                    f"{row.species}\t{row.alg}\t{row.best_chromosome}\t"
                    f"{row.alg_gene_count}\t{row.chromosome_gene_count}\t"
                    f"{row.alg_genes_on_best_chromosome}\t"
                    f"{row.splitting_parameter:.6f}\t"
                    f"{row.combining_parameter:.6f}\t"
                    f"{1 - row.splitting_parameter:.6f}\t"
                    f"{1 - row.combining_parameter:.6f}\t"
                    f"{row.rearrangement_index:.6f}\n"
                )
