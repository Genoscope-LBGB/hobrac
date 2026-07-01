"""Tile the per-reference ALG dotplots into a single high-resolution grid.

Each cell embeds one reference's finished dotplot PNG (axes already baked in by
dotplotrs) at native resolution. Every cell carries a pale-blue header strip with
the reference (X-axis) name, and each grid row carries a pale-blue vertical panel
on its left with the assembly (Y-axis) name. The grid lets the user see every
reference at once while staying zoomable for exploration.

Reference title resolution (first hit wins):
  1. --jcvi-names, when a full ordered list is given (1 assembly + N references,
     matching the karyotype labels). A single name applies only to the assembly,
     so references fall through.
  2. ``# Organism name:`` from reference/<accession>_assembly_report.txt.
  3. the accession itself.

Assembly title: the first --jcvi-names entry when any names are given, else the
``--assembly-name`` value (hobrac's ``-n/--name``).
"""

import argparse
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

# resolve_titles / resolve_assembly_title live in io.py (no matplotlib) and are
# re-exported here for the grid CLI and the karyotype labels in output.py.
from .io import (  # noqa: E402,F401
    parse_organism_name,
    resolve_assembly_title,
    resolve_titles,
)


def read_accession_order(path):
    """Read accessions (one per line) in ranking order."""
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def grid_dims(n):
    """Square-ish grid: cols = ceil(sqrt(n)), rows = ceil(n / cols)."""
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)
    return rows, cols


# Header (reference) strip height and left (assembly) panel width, as fractions
# of a single cell's native image size.
HEADER_FRAC = 0.09
LEFT_FRAC = 0.05

# Panel text height as a fraction of the header strip, when auto-sizing the font.
FONT_FRAC = 0.55

# Gap (as a fraction of a cell's width) between the left assembly panel and the
# dotplots, so the panel doesn't sit flush against the plots.
LEFT_PAD_FRAC = 0.02

# Panel fill + text colours per theme.
PANEL_COLORS = {
    "light": {"panel": "#d6e6f4", "text": "#1a1a1a", "bg": "white"},
    "dark": {"panel": "#2b3a55", "text": "#e8eef7", "bg": "#111111"},
}

# Margins (px) baked into each dotplot PNG by dotplotrs, i.e. the empty border
# between the image edge and the plotted area. Used to line the header/left
# panels up with the plot box rather than the full image.
DOTPLOT_MARGINS = {"top": 40, "bottom": 170, "left": 120, "right": 10}


def render_grid(
    images,
    titles,
    assembly_title,
    output_path,
    theme="light",
    dpi=100,
    title_fontsize=None,
    margins=None,
):
    """Tile *images* into one PNG sized to native resolution.

    Each cell gets a pale-blue header strip carrying its reference *title*; each
    grid row gets a pale-blue vertical panel on the left carrying *assembly_title*.
    Both panels are sized to the dotplot's plotted area (inside *margins*) so they
    line up with the plot box rather than the full image.
    """
    margins = margins or DOTPLOT_MARGINS
    n = len(images)
    rows, cols = grid_dims(n)
    colors = PANEL_COLORS.get(theme, PANEL_COLORS["light"])

    # Size each cell to the largest source image so no cell is downsampled.
    cell_h = max(img.shape[0] for img in images)
    cell_w = max(img.shape[1] for img in images)
    header_h = max(1, round(cell_h * HEADER_FRAC))
    left_w = max(1, round(cell_w * LEFT_FRAC))

    # Auto-size the panel font to a readable fraction of the header strip (in
    # points: px * 72 / dpi). Widen the narrow left panel if needed so the
    # rotated assembly text still fits inside its blue strip.
    if title_fontsize is None:
        fs_px = round(header_h * FONT_FRAC)
        title_fontsize = fs_px * 72 / dpi
        left_w = max(left_w, round(fs_px * 1.6))

    # Padding between the left panel and the dotplots; content (headers + images)
    # starts at content_x, so headers and images shift together and stay aligned.
    left_pad = round(cell_w * LEFT_PAD_FRAC)
    content_x = left_w + left_pad

    block_h = header_h + cell_h  # one grid row = header strip + dotplot
    total_w = content_x + cols * cell_w
    total_h = rows * block_h

    fig = plt.figure(figsize=(total_w / dpi, total_h / dpi), dpi=dpi)
    fig.patch.set_facecolor(colors["bg"])

    def rect(x_px, top_px, w_px, h_px):
        """A pixel box measured from the top-left -> figure-fraction rect."""
        return [
            x_px / total_w,
            1 - (top_px + h_px) / total_h,
            w_px / total_w,
            h_px / total_h,
        ]

    def panel(box, text, rotation=0):
        ax = fig.add_axes(box)
        ax.set_facecolor(colors["panel"])
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_edgecolor(colors["panel"])
        ax.text(
            0.5,
            0.5,
            text,
            ha="center",
            va="center",
            rotation=rotation,
            fontsize=title_fontsize,
            color=colors["text"],
            transform=ax.transAxes,
        )

    # Plot-box extents inside a cell (offsets/size relative to the cell image).
    plot_w = cell_w - margins["left"] - margins["right"]
    plot_h = cell_h - margins["top"] - margins["bottom"]

    for i, (img, title) in enumerate(zip(images, titles)):
        r, c = divmod(i, cols)
        block_top = r * block_h
        x_left = content_x + c * cell_w

        # Header strip aligned to the plot box's x-range (over the plotted area).
        panel(rect(x_left + margins["left"], block_top, plot_w, header_h), title)

        iax = fig.add_axes(rect(x_left, block_top + header_h, cell_w, cell_h))
        iax.imshow(img, interpolation="none", aspect="auto")
        iax.axis("off")

    # Assembly panel: one per grid row, down the far left, aligned to the plot
    # box's y-range (the plotted area, inside the top/bottom margins).
    for r in range(rows):
        plot_top = r * block_h + header_h + margins["top"]
        panel(rect(0, plot_top, left_w, plot_h), assembly_title, rotation=90)

    fig.savefig(output_path, dpi=dpi, facecolor=fig.get_facecolor())
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
    parser.add_argument("--assembly-name", default="assembly", help="Assembly name for the left panels (hobrac -n/--name)")
    parser.add_argument("-o", "--output", required=True, help="Output grid PNG path")
    args = parser.parse_args()

    suffix = "_dark.png" if args.theme == "dark" else ".png"

    accessions = read_accession_order(args.accession_order)
    titles_all = resolve_titles(accessions, args.jcvi_names, args.reference_dir)
    assembly_title = resolve_assembly_title(args.jcvi_names, args.assembly_name)

    images, titles = [], []
    for accession, title in zip(accessions, titles_all):
        png = os.path.join(args.dotplots_dir, f"{accession}{suffix}")
        if not os.path.isfile(png):
            continue
        images.append(mpimg.imread(png))
        titles.append(title)

    if not images:
        raise SystemExit(f"No dotplot PNGs found in {args.dotplots_dir} for theme {args.theme}")

    render_grid(images, titles, assembly_title, args.output, theme=args.theme)


if __name__ == "__main__":
    main()
