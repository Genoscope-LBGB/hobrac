"""Best-effort chromosome renaming for manually provided reference genomes.

When HoBRAC downloads a reference itself, ``find_reference_genomes`` renames the
sequences to ``chr<name>`` using the NCBI assembly report. References provided
manually with ``-r`` never go through that step, so their original (often opaque)
sequence ids are kept. As a best effort, we scan each FASTA header for a
``chr<token>`` pattern (e.g. ``chr1``, ``chrX``, ``chr2L``, ``chrMT``) and, when
found, rename the sequence to that token. A per-reference TSV mapping old ids to
new ids is written so the renaming stays traceable.
"""

import re
import sys

# Match a ``chr`` token such as chr1, chrX, chr2L, chrMT, chrUn... anywhere in a
# header. The character right after ``chr`` must be a digit or one of the usual
# chromosome letters; this deliberately excludes the word "chromosome" (chr + o).
# The token must not be glued to a preceding word character so we pick up whole
# tokens only.
CHR_PATTERN = re.compile(r"(?<![A-Za-z0-9_])chr[0-9XYZWMU][A-Za-z0-9_]*", re.IGNORECASE)


def find_chr_name(header: str):
    """Return the first ``chr<token>`` found in a FASTA header, or None."""
    match = CHR_PATTERN.search(header)
    return match.group(0) if match else None


def rename_reference(src_path: str, dest_fasta: str, mapping_path: str):
    """Copy ``src_path`` to ``dest_fasta`` renaming chromosome sequences.

    For every sequence whose header contains a ``chr<token>`` pattern, the
    sequence is renamed to that token. Sequences without a match, or whose target
    name would collide with one already used, keep their original id. A TSV
    mapping (``old_name<TAB>new_name``) is written to ``mapping_path`` with one
    row per sequence.
    """
    mapping = []
    used = set()

    with open(src_path) as src, open(dest_fasta, "w") as dst:
        for line in src:
            if not line.startswith(">"):
                dst.write(line)
                continue

            header = line[1:].rstrip("\n")
            old_name = header.split()[0] if header.split() else ""
            candidate = find_chr_name(header)

            if candidate and candidate != old_name and candidate in used:
                print(
                    f"Warning: cannot rename '{old_name}' to '{candidate}' in"
                    f" {src_path}: name already used. Keeping original.",
                    file=sys.stderr,
                )
                candidate = None

            new_name = candidate if candidate else old_name
            used.add(new_name)
            mapping.append((old_name, new_name))

            if new_name != old_name:
                dst.write(f">{new_name}\n")
            else:
                dst.write(line)

    with open(mapping_path, "w") as out:
        print("old_name\tnew_name", file=out)
        for old_name, new_name in mapping:
            print(f"{old_name}\t{new_name}", file=out)

    return mapping
