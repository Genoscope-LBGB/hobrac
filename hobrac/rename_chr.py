"""Best-effort chromosome renaming for manually provided reference genomes.

When HoBRAC downloads a reference itself, ``find_reference_genomes`` renames the
sequences to ``chr<name>`` using the NCBI assembly report. References provided
manually with ``-r`` never go through that step, so their original (often opaque)
sequence ids are kept. As a best effort, we scan each FASTA header for a
``chr<token>`` pattern (e.g. ``chr1``, ``chrX``, ``chr2L``, ``chrMT``) and, when
found, rename the sequence to that token. A per-reference TSV mapping old ids to
new ids is written so the renaming stays traceable.
"""

import os
import re
import sys

from xopen import xopen

# A chromosome-like token: a number (optionally with a trailing arm letter such
# as 2L), a single sex/special letter (X, Y, Z, W, U), or MT/Un. This excludes
# assembly names like "GRCh38" so the Ensembl "chromosome:GRCh38:1" form is not
# misread.
CHR_TOKEN = r"(?:[0-9]+[A-Za-z]?|MT|Un|[XYZWU])"

# Match a literal ``chr`` token such as chr1, chrX, chr2L, chrMT, chrUn...
# anywhere in a header. The token must not be glued to a preceding word
# character so we pick up whole tokens only.
CHR_PATTERN = re.compile(
    r"(?<![A-Za-z0-9_])chr" + CHR_TOKEN + r"[A-Za-z0-9_]*", re.IGNORECASE
)

# Match the descriptive form "chromosome 1", "chromosome: 1", "chromosome 2L"...
# A single separator (colon or whitespace) must follow "chromosome" so plurals
# ("chromosomes") are ignored, and the token must be chromosome-like so the
# Ensembl "chromosome:GRCh38:1" form (assembly name right after the colon) does
# not match.
CHROMOSOME_PATTERN = re.compile(
    r"chromosome[:\s]\s*(" + CHR_TOKEN + r")", re.IGNORECASE
)


def fasta_basename(path: str) -> str:
    """Basename of a FASTA path without its extension.

    A trailing ``.gz`` is stripped first so ``ref.fasta.gz`` and ``ref.fasta``
    both yield ``ref``. Used to derive stable reference ids; main.py and the
    select_references checkpoint must agree, so both call this.
    """
    name = os.path.basename(path)
    if name.endswith(".gz"):
        name = name[:-3]
    return os.path.splitext(name)[0]


def find_chr_name(header: str):
    """Return a ``chr<token>`` name derived from a FASTA header, or None.

    A literal ``chr<token>`` already present in the header is returned as-is.
    Otherwise a descriptive ``chromosome <token>`` is normalized to
    ``chr<token>``.
    """
    match = CHR_PATTERN.search(header)
    if match:
        return match.group(0)

    match = CHROMOSOME_PATTERN.search(header)
    if match:
        return "chr" + match.group(1)

    return None


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

    with xopen(src_path) as src, open(dest_fasta, "w") as dst:
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
