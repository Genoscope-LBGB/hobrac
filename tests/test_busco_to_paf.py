"""Tests for busco_to_paf PAF generation, including the co:Z: color tag."""

from hobrac.busco_to_paf import generate_paf


def test_generate_paf_appends_co_tag(tmp_path):
    busco1 = {  # assembly (query)
        "123at4751": {"chr": "chr1", "start": 100, "end": 200},
        "9at4751": {"chr": "chr2", "start": 50, "end": 90},
        "absent": {"chr": "chr3", "start": 0, "end": 10},
    }
    busco2 = {  # reference (target)
        "123at4751": {"chr": "scfA", "start": 300, "end": 400},
        "9at4751": {"chr": "scfB", "start": 500, "end": 540},
    }
    len_query = {"chr1": 1000, "chr2": 800, "chr3": 600}
    len_target = {"scfA": 1200, "scfB": 900}

    out = tmp_path / "aln_busco.paf"
    generate_paf(busco1, busco2, len_query, len_target, str(out))

    lines = [line for line in out.read_text().splitlines() if line]
    # Only the two shared genes are written ("absent" is dropped).
    assert len(lines) == 2

    for line in lines:
        fields = line.split("\t")
        # 12 mandatory PAF fields + 1 optional tag
        assert len(fields) == 13
        tag = fields[12]
        assert tag.startswith("co:Z:")
        busco_id = tag[len("co:Z:"):]
        assert busco_id in busco2  # tag carries the shared gene id

    # Spot-check the full first line for the 123at4751 gene.
    first = next(line for line in lines if "co:Z:123at4751" in line)
    f = first.split("\t")
    assert f[0] == "chr1"      # query chr (assembly)
    assert f[5] == "scfA"      # target chr (reference)
    assert f[12] == "co:Z:123at4751"
