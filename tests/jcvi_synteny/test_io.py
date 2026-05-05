from hobrac.jcvi_synteny.io import parse_custom_colors


def test_parse_custom_colors_accepts_rgb_and_hex(tmp_path):
    color_file = tmp_path / "colors.tsv"
    color_file.write_text(
        "# BUSCO_ID color ALG\n"
        "g1\t141,78,106\tALG_A\n"
        "g2\t#00FF7F\tALG_B\n"
        "g3\t336699\tALG_C\n"
    )

    assert parse_custom_colors(str(color_file)) == {
        "g1": "#8d4e6a",
        "g2": "#00ff7f",
        "g3": "#336699",
    }


def test_parse_custom_colors_skips_malformed_colors(tmp_path):
    color_file = tmp_path / "colors.tsv"
    color_file.write_text(
        "good\t0,128,255\tALG_A\n"
        "too_high\t256,0,0\tALG_B\n"
        "too_few\t1,2\tALG_C\n"
        "bad_hex\t#12xx56\tALG_D\n"
        "short_hex\t#abc\tALG_E\n"
    )

    assert parse_custom_colors(str(color_file)) == {"good": "#0080ff"}
