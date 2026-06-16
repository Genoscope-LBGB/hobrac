"""Tile the per-reference ALG dotplots into a single high-resolution grid.

Each cell embeds one reference's finished dotplot PNG (axes already baked in by
dotplotrs) at native resolution, with a species title on top. The grid lets the
user see every reference at once while staying zoomable for exploration.

Title resolution per reference (first hit wins):
  1. --jcvi-names, when a full ordered list is given (1 assembly + N references,
     matching the karyotype labels). A single name applies only to the assembly,
     so references fall through.
  2. ``# Organism name:`` from reference/<accession>_assembly_report.txt.
  3. the accession itself.
"""

import argparse
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


def read_accession_order(path):
    """Read accessions (one per line) in ranking order."""
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def parse_organism_name(report_path):
    """Return the ``# Organism name:`` value from an NCBI assembly report.

    Strips the trailing common-name parenthetical (e.g. "Homo sapiens (human)"
    -> "Homo sapiens"). Returns "" if the file is missing or has no such line.
    """
    if not os.path.isfile(report_path):
        return ""
    with open(report_path) as f:
        for line in f:
            if line.startswith("# Organism name:"):
                name = line.split(":", 1)[1].strip()
                if "(" in name:
                    name = name.split("(", 1)[0].strip()
                return name
    return ""


def resolve_titles(accessions, jcvi_names, reference_dir):
    """Resolve one display title per reference accession (see module docstring)."""
    names_list = [n.strip() for n in jcvi_names.split(",")] if jcvi_names else []
    # A full list is 1 assembly + N references; names_list[1:] maps to references
    # in accession order. Anything else (empty / assembly-only / mismatched) falls
    # through to the report/accession tiers.
    ref_names = (
        names_list[1:] if len(names_list) == len(accessions) + 1 else [""] * len(accessions)
    )

    titles = []
    for accession, name in zip(accessions, ref_names):
        if name:
            titles.append(name)
            continue
        report = os.path.join(reference_dir, f"{accession}_assembly_report.txt")
        titles.append(parse_organism_name(report) or accession)
    return titles


def grid_dims(n):
    """Square-ish grid: cols = ceil(sqrt(n)), rows = ceil(n / cols)."""
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)
    return rows, cols


def render_grid(images, titles, output_path, dpi=100, title_fontsize=14):
    """Tile *images* with *titles* into one PNG sized to native resolution."""
    n = len(images)
    rows, cols = grid_dims(n)

    # Size each cell to the largest source image so no cell is downsampled.
    cell_h = max(img.shape[0] for img in images)
    cell_w = max(img.shape[1] for img in images)

    fig, axes = plt.subplots(
        rows,
        cols,
        figsize=(cols * cell_w / dpi, rows * cell_h / dpi),
        dpi=dpi,
    )
    axes = list(axes.flat) if n > 1 else [axes]

    for ax, img, title in zip(axes, images, titles):
        ax.imshow(img, interpolation="none")
        ax.set_title(title, fontsize=title_fontsize)
        ax.axis("off")
    for ax in axes[n:]:
        ax.axis("off")

    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Tile per-reference ALG dotplots into a titled grid."
    )
    parser.add_argument("--dotplots-dir", required=True, help="Directory of per-reference dotplot PNGs")
    parser.add_argument("--accession-order", required=True, help="Accessions in ranking order, one per line")
    parser.add_argument("--theme", choices=["light", "dark"], default="light")
    parser.add_argument("--reference-dir", default="reference", help="Directory holding *_assembly_report.txt files")
    parser.add_argument("--jcvi-names", default="", help="Comma-separated track names (1 assembly + N references)")
    parser.add_argument("-o", "--output", required=True, help="Output grid PNG path")
    args = parser.parse_args()

    suffix = "_dark.png" if args.theme == "dark" else ".png"

    accessions = read_accession_order(args.accession_order)
    titles_all = resolve_titles(accessions, args.jcvi_names, args.reference_dir)

    images, titles = [], []
    for accession, title in zip(accessions, titles_all):
        png = os.path.join(args.dotplots_dir, f"{accession}{suffix}")
        if not os.path.isfile(png):
            continue
        images.append(mpimg.imread(png))
        titles.append(title)

    if not images:
        raise SystemExit(f"No dotplot PNGs found in {args.dotplots_dir} for theme {args.theme}")

    render_grid(images, titles, args.output)


if __name__ == "__main__":
    main()
