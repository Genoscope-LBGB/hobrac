"""Regression coverage for karyotype labels and unambiguous ALG legends."""

import matplotlib.pyplot as plt
import pytest

from hobrac.jcvi_synteny import legend


def test_draw_track_labels_without_floor_keeps_base_size():
    fig = plt.figure(figsize=(10, 10), dpi=100)
    try:
        ax = fig.add_axes([0, 0, 1, 1])
        fig.canvas.draw()

        legend.draw_track_labels(
            ax,
            ["Short"],
            [(0.5, 0.1)],
            1000,
            1000,
            100,
            fig.canvas.get_renderer(),
            min_species_label_font_size=None,
        )

        assert len(ax.texts) == 1
        assert ax.texts[0].get_fontsize() == pytest.approx(12.96)
    finally:
        plt.close(fig)


def test_draw_track_labels_wrappable_overflow_honors_floor():
    fig = plt.figure(figsize=(10, 10), dpi=100)
    try:
        ax = fig.add_axes([0, 0, 1, 1])
        fig.canvas.draw()

        legend.draw_track_labels(
            ax,
            [
                "Extremely long species name that cannot fit "
                "in the narrow left label margin"
            ],
            [(0.5, 0.1)],
            1000,
            1000,
            100,
            fig.canvas.get_renderer(),
            min_species_label_font_size=20.0,
        )

        assert len(ax.texts) == 1
        text = ax.texts[0]
        assert text.get_fontsize() == 20.0
        assert "\n" in text.get_text()
        line_widths = [
            fig.canvas.get_renderer().get_text_width_height_descent(
                line, text.get_fontproperties(), False
            )[0]
            for line in text.get_text().splitlines()
        ]
        assert all(width > 92 for width in line_widths)
    finally:
        plt.close(fig)


def test_draw_track_labels_single_token_overflow_honors_floor():
    fig = plt.figure(figsize=(10, 10), dpi=100)
    try:
        ax = fig.add_axes([0, 0, 1, 1])
        fig.canvas.draw()

        legend.draw_track_labels(
            ax,
            ["X" * 100],
            [(0.5, 0.1)],
            1000,
            1000,
            100,
            fig.canvas.get_renderer(),
            min_species_label_font_size=20.0,
        )

        assert len(ax.texts) == 1
        text = ax.texts[0]
        assert text.get_fontsize() == 20.0
        assert "\n" not in text.get_text()
    finally:
        plt.close(fig)


def test_collect_legend_uses_single_display_chain_for_ambiguous_gene(tmp_path):
    gene_chains = tmp_path / "gene_chains.tsv"
    gene_chains.write_text(
        "chain_ids\tcustom_alg_id\tgene\tcolor\tdisplay_chain_id\n"
        "0,2\t-\tg1\t#123456\t2\n"
        "2\t-\tg2\t#123456\t2\n"
        "0\t-\tg3\tlightgrey\t0\n"
    )

    assert legend.collect_legend_entries(gene_chains) == [("2", "#123456")]


def test_collect_legend_preserves_custom_labels_with_or_without_chains(tmp_path):
    gene_chains = tmp_path / "gene_chains.tsv"
    gene_chains.write_text(
        "chain_ids\tcustom_alg_id\tgene\tcolor\tdisplay_chain_id\n"
        "1\tO2\tg1\t#abcdef\t1\n"
        "-\tO10\tg2\t#fedcba\t-\n"
        "-\t-\tg3\tlightgrey\t-\n"
    )

    assert legend.collect_legend_entries(gene_chains) == [
        ("O2", "#abcdef"),
        ("O10", "#fedcba"),
    ]
