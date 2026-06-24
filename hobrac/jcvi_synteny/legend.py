"""Composite track labels and an ALG colour legend onto the karyotype PNG.

The karyotype image is rendered upstream by ``jcvi.graphics.karyotype`` with its
own label column left blank (long species names would otherwise overlap the
tracks); this step reads the finished ``karyotype.png`` and draws the species
labels right-aligned in the left margin (from the ``# track_labels`` comment in
the layouts file), plus a legend panel of stacked colour bricks on the right, one
brick per ALG actually shown in the plot.

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

Bricks are stacked top-to-bottom in natural-sorted label order (``A1``, ``A2``,
``B1`` ... / numeric chain ids), confined to the vertical band of the plotted
tracks (top track's labels down to the bottom track's chromosomes, read from the
``layouts`` file) rather than the full image height. If nothing is coloured, the
karyotype is left untouched and no panel is drawn.
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
PANEL_WIDTH_FRAC = 0.10
MAX_BRICK_FRAC = 0.055

# Brick label height as a fraction of the brick height, and the title strip
# height as a fraction of the image height.
LABEL_FRAC = 0.42
TITLE_FRAC = 0.035

# Padding (fraction of image width) between the karyotype and the bricks, and
# inside the panel on either side of a brick.
GAP_FRAC = 0.012
BRICK_INSET_FRAC = 0.18

# The legend is confined to the vertical band of the actual tracks: from a little
# above the top track's bar (to reach its labels) down to the bottom track's bar
# (its chromosomes). These pads, in fractions of image height, set how far above
# the top bar and below the bottom bar the band reaches.
TOP_BAND_PAD_FRAC = 0.04
BOTTOM_BAND_PAD_FRAC = 0.0

# Track labels (drawn here instead of by jcvi so long names cannot overlap the
# tracks): base font height as a fraction of the image height (so it does not
# balloon when there are few tracks), capped at a fraction of the inter-track
# spacing so two wrapped lines still fit between adjacent tracks. LABEL_GAP_FRAC
# is the gap (fraction of image width) between a label's right edge and its track.
LABEL_FONT_FRAC = 0.018
LABEL_MAX_SPACING_FRAC = 0.34
LABEL_GAP_FRAC = 0.008


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


def read_track_band(layouts_path):
    """Return ``(band_top, band_bottom)`` as fractions of image height (top=0).

    Reads the karyotype ``layouts`` file, whose track rows start with the track's
    ``y`` coordinate in jcvi axes space (1 = top of figure). The top track has the
    largest ``y`` and the bottom track the smallest; jcvi draws a track's bar at
    pixel fraction ``1 - y`` from the top. The band runs from a little above the
    top bar (to include the top track's labels) down to the bottom bar. Returns
    ``None`` if no track rows are found (caller falls back to the full height).
    """
    ys = []
    with open(layouts_path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            first = stripped.split(",")[0].strip()
            try:
                ys.append(float(first))
            except ValueError:
                continue  # edge rows ("e, 0, 1, links.simple") have no leading y
    if not ys:
        return None
    top_y, bottom_y = max(ys), min(ys)
    band_top = max(0.0, (1.0 - top_y) - TOP_BAND_PAD_FRAC)
    band_bottom = min(1.0, (1.0 - bottom_y) + BOTTOM_BAND_PAD_FRAC)
    return band_top, band_bottom


def read_track_labels(layouts_path):
    """Return the per-track display labels from the ``# track_labels`` comment.

    ``generate_layouts_file`` stashes the labels there (tab-separated, in track
    order) because jcvi's own label column is left blank. Returns ``[]`` when the
    comment is absent.
    """
    with open(layouts_path) as f:
        for line in f:
            if line.startswith("# track_labels"):
                return line.rstrip("\n").split("\t")[1:]
    return []


def read_track_positions(layouts_path):
    """Return ``[(y, xstart)]`` per track row, in file order.

    ``y`` is the jcvi axes y (1 = top of figure) and ``xstart`` the track's left
    edge (fraction of image width). Comment and edge rows are skipped.
    """
    rows = []
    with open(layouts_path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            cols = [c.strip() for c in stripped.split(",")]
            try:
                rows.append((float(cols[0]), float(cols[1])))
            except (ValueError, IndexError):
                continue  # edge rows ("e, 0, 1, links.simple")
    return rows


def _wrap_two_lines(text):
    """Split a multi-word *text* into two balanced lines on whitespace.

    Returns ``None`` for a single token (e.g. an accession ``CA_96446885``): a
    mid-token split reads worse than shrinking it onto one line, which the caller
    does instead.
    """
    words = text.split()
    if len(words) < 2:
        return None
    k = min(
        range(1, len(words)),
        key=lambda i: abs(len(" ".join(words[:i])) - len(" ".join(words[i:]))),
    )
    return " ".join(words[:k]) + "\n" + " ".join(words[k:])


def draw_track_labels(ax, labels, positions, img_w, img_h, dpi, renderer):
    """Draw right-aligned track labels in the left margin of the karyotype.

    *ax* spans the karyotype image with ``xlim(0, 1)`` / ``ylim(1, 0)`` so jcvi's
    coordinates map straight through: a track at jcvi-y ``y`` has its bar at axes
    ``1 - y`` and starts at axes-x ``xstart``. Each label is right-aligned just
    left of its track; if a single line is wider than the available margin it is
    wrapped to two lines, then shrunk to fit. ponytail: font auto-shrink guarantees
    the label always fits rather than running off the image edge.
    """
    n = len(positions)
    spacing = (0.8 / max(n - 1, 1)) * img_h  # px between adjacent tracks
    base_fs_px = min(LABEL_FONT_FRAC * img_h, LABEL_MAX_SPACING_FRAC * spacing)
    base_fs = base_fs_px * 72 / dpi  # px -> points

    def width_pt(text, fs):
        t = ax.text(0, 0, text, fontsize=fs)
        w = t.get_window_extent(renderer).width
        t.remove()
        return w

    for (y, xstart), text in zip(positions, labels):
        if not text:
            continue
        x_right = xstart - LABEL_GAP_FRAC
        max_w = x_right * img_w  # display px available to the left edge
        rendered, fs = text, base_fs
        if width_pt(rendered, fs) > max_w:
            wrapped = _wrap_two_lines(text)  # None for single-token names
            if wrapped is not None:
                rendered = wrapped
            w = width_pt(rendered, fs)
            if w > max_w:
                fs *= max_w / w  # shrink (wrapped or single line) to fit the margin
        ax.text(
            x_right,
            1 - y,
            rendered,
            ha="right",
            va="center",
            fontsize=fs,
            color="#000000",
            linespacing=0.95,
        )


def _text_color(face_color):
    """Black or white, whichever reads better on *face_color* (sRGB luminance)."""
    r, g, b = mcolors.to_rgb(face_color)
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    return "#000000" if luminance > 0.55 else "#ffffff"


def render_legend(
    karyotype_path,
    entries,
    output_path,
    title="ALGs",
    dpi=100,
    band=(0.0, 1.0),
    labels=None,
    positions=None,
):
    """Composite track labels (left) and an ALG brick legend (right) onto the PNG.

    Track *labels* are drawn right-aligned in the karyotype's left margin at the
    positions read from the layouts file (see ``draw_track_labels``). The ALG
    legend is only added when *entries* is non-empty; *band* is ``(top, bottom)``
    in fractions of image height (top = 0), confining the bricks to the plotted
    tracks. With no entries the figure keeps the karyotype's own width.
    """
    base = mpimg.imread(karyotype_path)
    img_h, img_w = base.shape[0], base.shape[1]

    has_panel = bool(entries)
    panel_w = max(1, round(img_w * PANEL_WIDTH_FRAC)) if has_panel else 0
    gap = round(img_w * GAP_FRAC) if has_panel else 0
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

    # Track labels overlay: an axes spanning the karyotype image, in jcvi's own
    # 0..1 coordinate space (y inverted) so labels land on the bars jcvi drew.
    if labels and positions:
        tax = fig.add_axes(rect(0, 0, img_w, img_h))
        tax.set_xlim(0, 1)
        tax.set_ylim(1, 0)
        tax.axis("off")
        draw_track_labels(
            tax, labels, positions, img_w, img_h, dpi, fig.canvas.get_renderer()
        )

    if has_panel:
        # Legend panel occupies the right column; one axes spanning it, drawn in
        # its own 0..1 coordinate space with y inverted so the first entry is on top.
        lax = fig.add_axes(rect(img_w + gap, 0, panel_w, img_h))
        lax.set_xlim(0, 1)
        lax.set_ylim(1, 0)
        lax.axis("off")

        n = len(entries)
        band_top, band_bottom = band
        title_h = TITLE_FRAC  # in axes-fraction units (1.0 == full image height)
        avail = (band_bottom - band_top) - title_h
        brick_h = min(MAX_BRICK_FRAC, avail / n)
        inset = BRICK_INSET_FRAC

        label_fs = round(brick_h * total_h * LABEL_FRAC) * 72 / dpi

        lax.text(
            0.5,
            band_top + title_h / 2,
            title,
            ha="center",
            va="center",
            fontweight="bold",
            fontsize=label_fs,
            color="#000000",
        )

        bricks_top = band_top + title_h
        for idx, (label, color) in enumerate(entries):
            top = bricks_top + idx * brick_h
            lax.add_patch(
                mpatches.Rectangle(
                    (inset, top),
                    1 - 2 * inset,
                    brick_h,
                    facecolor=color,
                    edgecolor="#000000",
                    linewidth=0.8,
                )
            )
            lax.text(
                0.5,
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
        "--layouts",
        default=None,
        help="jcvi layouts file; supplies the track labels (drawn in the left "
        "margin) and confines the ALG legend to the track band",
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

    band = (0.0, 1.0)
    labels, positions = [], []
    if args.layouts:
        track_band = read_track_band(args.layouts)
        if track_band:
            band = track_band
        labels = read_track_labels(args.layouts)
        positions = read_track_positions(args.layouts)

    if not entries and not (labels and positions):
        print("Nothing to composite; leaving karyotype unchanged.")
        return

    output_path = args.output or args.karyotype
    render_legend(
        args.karyotype,
        entries,
        output_path,
        title=args.title,
        band=band,
        labels=labels,
        positions=positions,
    )


if __name__ == "__main__":
    main()
