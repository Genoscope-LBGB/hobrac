"""Regression coverage for karyotype species-label sizing and CLI plumbing."""

import sys
from unittest.mock import patch

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


def test_render_legend_forwards_species_label_floor(tmp_path):
    numpy = pytest.importorskip("numpy")
    import matplotlib.image as mpimg

    input_png = tmp_path / "input.png"
    output_png = tmp_path / "output.png"
    mpimg.imsave(input_png, numpy.ones((100, 100, 3)))

    with patch.object(legend, "draw_track_labels") as draw_track_labels:
        legend.render_legend(
            input_png,
            [],
            output_png,
            labels=["Species One"],
            positions=[(0.5, 0.1)],
            min_species_label_font_size=14.5,
        )

    draw_track_labels.assert_called_once()
    _, kwargs = draw_track_labels.call_args
    assert kwargs == {"min_species_label_font_size": 14.5}
    assert output_png.is_file()


def test_main_forwards_species_label_floor(tmp_path, monkeypatch):
    gene_chains = tmp_path / "gene_chains.tsv"
    gene_chains.write_text("chain_id\tcustom_alg_id\tgene\tcolor\n")
    layouts = tmp_path / "layouts"
    layouts.write_text(
        "# track_labels\tSpecies One\n"
        "0.5, 0.1, 0.96, 0, black, , bottom, species.bed\n"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "karyotype_legend",
            "--gene-chains",
            str(gene_chains),
            "--karyotype",
            str(tmp_path / "karyotype.png"),
            "--layouts",
            str(layouts),
            "--min-species-label-font-size",
            "14.5",
        ],
    )

    with patch.object(legend, "render_legend") as render_legend:
        legend.main()

    render_legend.assert_called_once()
    _, kwargs = render_legend.call_args
    assert kwargs["min_species_label_font_size"] == 14.5


@pytest.mark.parametrize("value", ["0", "-1", "nan"])
def test_main_rejects_non_positive_or_non_finite_label_floor(
    value, monkeypatch, capsys
):
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "karyotype_legend",
            "--gene-chains",
            "gene_chains.tsv",
            "--karyotype",
            "karyotype.png",
            "--min-species-label-font-size",
            value,
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        legend.main()

    assert exc_info.value.code == 2
    assert "must be a finite positive number" in capsys.readouterr().err
