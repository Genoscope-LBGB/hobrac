import argparse
import sys
from collections.abc import Iterable


def deduplicate(
    rows: Iterable[tuple[str, str, str, str]],
) -> list[tuple[str, str, str]]:
    """Keep one (accession, name, taxid) per assembly, in first-seen order.

    Rows are (Assembly Name, Accession, Paired Accession, Taxid). NCBI lists an
    assembly once per GenBank/RefSeq record and may keep an older GCA version
    when it is the one paired with the current GCF, so an assembly is identified
    by its unversioned GCA accession (its own or the paired one). Assembly names
    are not unique across assemblies and cannot serve as the key. The newest GCA
    version wins; a GCF row is kept only when the assembly has no GCA row.
    """
    best: dict[str, tuple[tuple[bool, int], str, str, str]] = {}
    for name, accession, paired, taxid in rows:
        is_gca = accession.startswith("GCA_")
        gca = accession if is_gca else paired
        key = (gca or accession).split(".")[0]
        rank = (is_gca, int(accession.rsplit(".", 1)[1]))
        if key not in best or rank > best[key][0]:
            best[key] = (rank, accession, name, taxid)
    return [(accession, name, taxid) for _, accession, name, taxid in best.values()]


def main():
    parser = argparse.ArgumentParser(
        prog="dedup_ncbi",
        description=(
            "Deduplicate an NCBI Datasets assembly TSV into the precompute_mash"
            " input format: GCA<TAB>Name<TAB>Taxid."
        ),
    )
    parser.add_argument("input", help="NCBI Datasets TSV (with header)")
    parser.add_argument(
        "-o", "--output", default="-", help="Output file (default: stdout)"
    )
    args = parser.parse_args()

    rows = []
    with open(args.input) as inf:
        next(inf, None)  # drop header
        for line in inf:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            rows.append(tuple(fields[:4]))

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    try:
        for accession, name, taxid in deduplicate(rows):
            print(f"{accession}\t{name}\t{taxid}", file=out)
    finally:
        if out is not sys.stdout:
            out.close()


if __name__ == "__main__":
    main()
