import warnings
from collections import defaultdict
from typing import Dict

from .models import BuscoGene


def read_fasta_sizes(fasta_path: str) -> Dict[str, int]:
    """
    Read sequence sizes from a fasta file.

    Args:
        fasta_path: Path to the fasta file

    Returns:
        Dictionary mapping sequence names to their lengths
    """
    sizes = {}
    current_name = None
    current_length = 0

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name is not None:
                    sizes[current_name] = current_length
                current_name = line[1:].split()[0]
                current_length = 0
            else:
                current_length += len(line)

        if current_name is not None:
            sizes[current_name] = current_length

    return sizes


def parse_custom_colors(color_file: str) -> Dict[str, str]:
    """
    Parse custom color file and normalize colors to hex.

    Color file format (tab or space separated):
    BUSCO_ID     COLOR        ALG_NAME
    10000at6447  141,78,106   A2
    10001at6447  #8d4e6a     A2
    10002at6447  8d4e6a      A2

    Args:
        color_file: Path to the custom color file

    Returns:
        Dictionary mapping BUSCO ID to hex color string
    """
    custom_colors = {}
    with open(color_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            busco_id = parts[0]
            color = _normalize_custom_color(parts[1])
            if color:
                custom_colors[busco_id] = color
    return custom_colors


def _normalize_custom_color(color: str) -> str:
    """Return a normalized #rrggbb color, or an empty string if invalid."""
    if "," in color:
        try:
            channels = [int(channel) for channel in color.split(",")]
        except ValueError:
            return ""
        if len(channels) != 3 or any(channel < 0 or channel > 255 for channel in channels):
            return ""
        return f"#{channels[0]:02x}{channels[1]:02x}{channels[2]:02x}"

    hex_color = color[1:] if color.startswith("#") else color
    if len(hex_color) == 6 and all(c in "0123456789abcdefABCDEF" for c in hex_color):
        return f"#{hex_color.lower()}"

    return ""


def parse_custom_algs(color_file: str) -> Dict[str, str]:
    """
    Parse ALG labels from a custom color file.

    Color file format (tab or space separated):
    BUSCO_ID     COLOR        ALG_NAME

    Args:
        color_file: Path to custom color file

    Returns:
        Dictionary mapping BUSCO ID to ALG label from the third column.
    """
    custom_algs = {}
    with open(color_file, "r") as f:
        for line_number, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 3:
                warnings.warn(
                    f"{color_file}:{line_number}: expected 3 columns "
                    f"(BUSCO_ID, COLOR, ALG_NAME), found {len(parts)}",
                    UserWarning,
                    stacklevel=2,
                )
            if len(parts) >= 3:
                custom_algs[parts[0]] = parts[2]
    return custom_algs


def read_busco_tsv(file_path: str, min_busco_genes: int = 0) -> Dict[str, BuscoGene]:
    """
    Parse BUSCO full_table.tsv and return only Complete single-copy genes.

    Args:
        file_path: Path to BUSCO full_table.tsv file
        min_busco_genes: Minimum complete BUSCO genes required per sequence.
                         Sequences with fewer genes are excluded.

    Returns:
        Dictionary mapping BUSCO ID to BuscoGene
    """
    busco_data = {}
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 5 and parts[1] == "Complete":
                busco_id = parts[0]
                busco_data[busco_id] = BuscoGene(
                    busco_id=busco_id,
                    chromosome=parts[2],
                    start=int(parts[3]),
                    end=int(parts[4]),
                )

    if min_busco_genes > 0:
        busco_data = filter_by_min_genes(busco_data, min_busco_genes)

    return busco_data


def filter_by_min_genes(
    busco_data: Dict[str, BuscoGene], min_genes: int
) -> Dict[str, BuscoGene]:
    """
    Filter BUSCO data to keep only sequences with at least min_genes.

    Args:
        busco_data: Dictionary mapping BUSCO ID to BuscoGene
        min_genes: Minimum number of genes required per sequence

    Returns:
        Filtered dictionary with only genes from qualifying sequences
    """
    seq_counts: Dict[str, int] = defaultdict(int)
    for gene in busco_data.values():
        seq_counts[gene.chromosome] += 1

    valid_seqs = {seq for seq, count in seq_counts.items() if count >= min_genes}

    return {
        busco_id: gene
        for busco_id, gene in busco_data.items()
        if gene.chromosome in valid_seqs
    }
