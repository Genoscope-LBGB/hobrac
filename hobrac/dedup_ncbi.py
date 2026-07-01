import argparse
import sys


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

    # Input cols: 1=Assembly Name, 2=Accession, 3=Paired Accession, 4=Taxid.
    # Each assembly appears twice (one GCA row, one GCF row) sharing the Assembly
    # Name. Keep one row per name, preferring the GCA accession; if a name only
    # has a GCF row, keep that.
    best = {}
    order = []
    with open(args.input) as inf:
        next(inf, None)  # drop header
        for line in inf:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            name, accession, _, taxid = fields[:4]
            is_gca = accession.startswith("GCA_")
            if name not in best:
                order.append(name)
                best[name] = (accession, taxid, is_gca)
            elif is_gca and not best[name][2]:
                best[name] = (accession, taxid, is_gca)

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    try:
        for name in order:
            accession, taxid, _ = best[name]
            print(f"{accession}\t{name}\t{taxid}", file=out)
    finally:
        if out is not sys.stdout:
            out.close()


if __name__ == "__main__":
    main()
