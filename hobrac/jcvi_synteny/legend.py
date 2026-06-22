"""Composite an ALG colour legend onto the right margin of the karyotype PNG.

The karyotype image is rendered upstream by ``jcvi.graphics.karyotype``; this
step reads the finished ``karyotype.png`` and draws a legend panel of stacked
colour bricks beside it, one brick per ALG actually shown in the plot.

Legend contents come straight from ``gene_chains.tsv`` (columns ``chain_id``,
``custom_alg_id``, ``gene``, ``color``):

  * Only genes drawn as coloured ribbons count: rows whose ``color`` is
    ``lightgrey`` (non-significant / unlisted) are ignored, so the legend lists
    only the ALGs the plot actually shows.
  * When a named colour scheme is active each coloured gene carries a
    ``custom_alg_id`` (e.g. ``O1``); the legend then shows one brick per ALG
    name. With default chain colouring ``custom_alg_id`` is ``-`` and the legend
    shows one brick per chain, labelled with the chain id -- even where the
    palette reuses a colour across chains, so the brick count stays faithful to
    the chain count.

Bricks are laid out in two columns, filled column-major in natural-sorted label
order (``A1``, ``A2``, ``B1`` ... / numeric chain ids) so entries read top-to-
bottom down the left column then down the right, and the whole block is centred
vertically on the page. If nothing is coloured, the karyotype is left untouched
and no panel is drawn.
"""

import argparse
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from .models import DEFAULT_COLOR  # noqa: E402

# Legend panel width as a fraction of the karyotype image width, and the largest
# a single brick is allowed to grow to (as a fraction of image height) so a
# short legend keeps tidy bricks instead of a few giant blocks.
PANEL_WIDTH_FRAC = 0.22
MAX_BRICK_FRAC = 0.055

# Brick label height as a fraction of the brick height, and the title strip
# height as a fraction of the image height.
LABEL_FRAC = 0.42
TITLE_FRAC = 0.035

# Padding (fraction of image width) between the karyotype and the bricks, and
# inside the panel on either side of a brick.
GAP_FRAC = 0.012
BRICK_INSET_FRAC = 0.18

# Gap between the two brick columns, as a fraction of the panel width.
COLUMN_GAP_FRAC = 0.10


def _natural_key(label):
    """Split digit runs so ``A2`` sorts before ``A10`` and ``2`` before ``10``."""
    return [
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r"(\d+)", label)
    ]


def collect_legend_entries(gene_chains_path):
    """Return ``[(label, hex_color)]`` for every ALG shown, in legend order.

    Reads ``gene_chains.tsv`` and keeps one entry per legend unit: per
    ``custom_alg_id`` when a named scheme is in use, else per ``chain_id``. Rows
    coloured ``lightgrey`` (the default, non-significant colour) are skipped. The
    first colour seen for a unit wins; entries come back natural-sorted by label.
    """
    entries = {}
    with open(gene_chains_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            i_chain = header.index("chain_id")
            i_alg = header.index("custom_alg_id")
            i_color = header.index("color")
        except ValueError as exc:
            raise SystemExit(
                f"{gene_chains_path}: missing expected column ({exc})"
            )
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= i_color:
                continue
            color = cols[i_color]
            if not color or color == DEFAULT_COLOR:
                continue
            alg = cols[i_alg]
            label = alg if alg and alg != "-" else cols[i_chain]
            if not label or label == "-":
                continue
            entries.setdefault(label, color)
    return sorted(entries.items(), key=lambda item: _natural_key(item[0]))


def _text_color(face_color):
    """Black or white, whichever reads better on *face_color* (sRGB luminance)."""
    r, g, b = mcolors.to_rgb(face_color)
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    return "#000000" if luminance > 0.55 else "#ffffff"


def render_legend(karyotype_path, entries, output_path, title="ALGs", dpi=100):
    """Composite *entries* as a brick legend onto the right of the karyotype PNG."""
    base = mpimg.imread(karyotype_path)
    img_h, img_w = base.shape[0], base.shape[1]

    panel_w = max(1, round(img_w * PANEL_WIDTH_FRAC))
    gap = round(img_w * GAP_FRAC)
    total_w = img_w + gap + panel_w
    total_h = img_h

    fig = plt.figure(figsize=(total_w / dpi, total_h / dpi), dpi=dpi)
    fig.patch.set_facecolor("white")

    def rect(x_px, top_px, w_px, h_px):
        """Pixel box measured from the top-left -> figure-fraction rect."""
        return [
            x_px / total_w,
            1 - (top_px + h_px) / total_h,
            w_px / total_w,
            h_px / total_h,
        ]

    # Karyotype on the left, at native pixel size so it is not resampled.
    kax = fig.add_axes(rect(0, 0, img_w, img_h))
    kax.imshow(base, interpolation="none", aspect="auto")
    kax.axis("off")

    # Legend panel occupies the right column; one axes spanning it, drawn in its
    # own 0..1 coordinate space with y inverted so the first entry sits on top.
    lax = fig.add_axes(rect(img_w + gap, 0, panel_w, img_h))
    lax.set_xlim(0, 1)
    lax.set_ylim(1, 0)
    lax.axis("off")

    n = len(entries)
    title_h = TITLE_FRAC  # in axes-fraction units (panel spans full height)
    inset = BRICK_INSET_FRAC

    # Two columns, filled column-major so natural-sorted entries read top-to-
    # bottom down the left column, then down the right. Bricks keep the size they
    # would have as one full-height column of *n* entries, so splitting into two
    # columns leaves the stack about half as tall -- which we then centre.
    rows_per_col = (n + 1) // 2
    avail = 1.0 - title_h
    brick_h = min(MAX_BRICK_FRAC, avail / n)

    # Two equal columns with a gap between them (x runs 0..1 across the panel).
    col_gap = COLUMN_GAP_FRAC
    col_w = (1.0 - col_gap) / 2
    col_x = [0.0, col_w + col_gap]  # left edge of each column

    # Centre the whole block (title + bricks) vertically in the panel.
    block_h = title_h + rows_per_col * brick_h
    top0 = (1.0 - block_h) / 2

    label_fs = round(brick_h * total_h * LABEL_FRAC) * 72 / dpi

    lax.text(
        0.5,
        top0 + title_h / 2,
        title,
        ha="center",
        va="center",
        fontweight="bold",
        fontsize=label_fs,
        color="#000000",
    )

    bricks_top = top0 + title_h
    for idx, (label, color) in enumerate(entries):
        col, row = divmod(idx, rows_per_col)
        x0 = col_x[col]
        top = bricks_top + row * brick_h
        lax.add_patch(
            mpatches.Rectangle(
                (x0 + inset * col_w, top),
                col_w * (1 - 2 * inset),
                brick_h,
                facecolor=color,
                edgecolor="#000000",
                linewidth=0.8,
            )
        )
        lax.text(
            x0 + col_w / 2,
            top + brick_h / 2,
            label,
            ha="center",
            va="center",
            fontsize=label_fs,
            color=_text_color(color),
        )

    fig.savefig(output_path, dpi=dpi, facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Composite an ALG colour legend onto the karyotype PNG."
    )
    parser.add_argument(
        "--gene-chains", required=True, help="Path to gene_chains.tsv"
    )
    parser.add_argument(
        "--karyotype", required=True, help="Karyotype PNG to annotate (overwritten)"
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output PNG (defaults to overwriting --karyotype)",
    )
    parser.add_argument("--title", default="ALGs", help="Legend header text")
    args = parser.parse_args()

    entries = collect_legend_entries(args.gene_chains)
    if not entries:
        print("No coloured ALGs in gene_chains.tsv; leaving karyotype unchanged.")
        return

    output_path = args.output or args.karyotype
    render_legend(args.karyotype, entries, output_path, title=args.title)


if __name__ == "__main__":
    main()
