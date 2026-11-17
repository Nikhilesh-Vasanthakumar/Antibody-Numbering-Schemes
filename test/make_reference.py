"""Generate test/fixtures/reference.json from the Python `antibody-numbering` package.

The Python package is cross-checked against ANARCI. This script turns its
alignments into the Stockholm-style inputs these JS functions take (RF line +
aligned row, with flanking residues shown as insert columns the way hmmalign
shows them) and stores the expected numbering for every scheme.

Usage (in an environment with antibody-numbering installed):
    python test/make_reference.py sequences.fasta[.gz] test/fixtures/reference.json [--diverse N]

--diverse N keeps at most N domains, one per distinct loop-length/insertion
pattern, so the fixture stays small but covers the edge cases.
"""

import argparse
import gzip
import json
from pathlib import Path

from antibody_numbering.alignment import find_domains
from antibody_numbering.schemes import NumberingError, apply_scheme

VALID = set("ACDEFGHIKLMNPQRSTVWY")
SCHEMES = ("imgt", "kabat", "chothia", "aho")


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    name, buf = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name:
        yield name, "".join(buf)


def stockholm_row(states, seq):
    """RF line and aligned row over all 128 match states, flanks as insert columns."""
    numbered_idx = [si for _, si in states if si is not None]
    first_si, last_si = numbered_idx[0], numbered_idx[-1]
    ann, row = [], []
    for aa in seq[:first_si]:  # N-terminal flank
        ann.append(".")
        row.append(aa.lower())
    current = 0
    for (state, kind), si in states:
        if kind == "i":
            ann.append(".")
            row.append(seq[si].lower())
            continue
        while current < state - 1:  # unaligned states before this one
            current += 1
            ann.append("x")
            row.append("-")
        current = state
        ann.append("x")
        row.append("-" if kind == "d" else seq[si])
    while current < 128:
        current += 1
        ann.append("x")
        row.append("-")
    for aa in seq[last_si + 1 :]:  # C-terminal flank
        ann.append(".")
        row.append(aa.lower())
    return "".join(ann), "".join(row)


def expected(states, seq, scheme, chain):
    try:
        numbered, _, _ = apply_scheme(states, seq, scheme, chain)
    except NumberingError as exc:
        return {"error": str(exc)}
    return {
        "concatSequence": "".join(aa for _, aa in numbered),
        "numbering": [n if not ins.strip() else f"{n}{ins.strip()}" for (n, ins), _ in numbered],
    }


def signature(rec):
    """Loop lengths and insertion pattern: what makes a domain a distinct test case."""
    parts = [rec["chain"]]
    for scheme in SCHEMES:
        exp = rec["expected"].get(scheme)
        if not exp:
            continue
        if "error" in exp:
            parts.append(f"{scheme}:err")
        else:
            parts.append(scheme + ":" + ",".join(str(p) for p in exp["numbering"] if isinstance(p, str)))
    imgt = rec["expected"]["imgt"]
    if "numbering" in imgt:
        nums = [int(str(p).rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")) for p in imgt["numbering"]]
        parts.append(f"start{nums[0]}end{nums[-1]}")
        for lo, hi in ((27, 38), (56, 65), (105, 117)):
            parts.append(str(sum(lo <= n <= hi for n in nums)))
    return "|".join(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--diverse", type=int, default=0)
    args = ap.parse_args()

    records = [(n, s) for n, s in read_fasta(args.fasta) if not set(s) - VALID]
    out = []
    for (name, seq), domains in zip(records, find_domains([s for _, s in records]), strict=True):
        if len(domains) != 1:
            continue  # the JS functions take one domain per row
        d = domains[0]
        ann, row = stockholm_row(d.states, seq)
        schemes = SCHEMES if d.chain_type in "HKL" else ("imgt", "aho")
        out.append({
            "name": name,
            "chain": d.chain_type,
            "annotations": ann,
            "sequence": row,
            "expected": {s: expected(d.states, seq, s, d.chain_type) for s in schemes},
        })
    if args.diverse:
        seen, kept = set(), []
        for rec in out:
            sig = signature(rec)
            if sig not in seen:
                seen.add(sig)
                kept.append(rec)
        out = kept[: args.diverse]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, separators=(",", ":")) + "\n", encoding="utf-8")
    print(f"wrote {len(out)} domains to {args.out}")


if __name__ == "__main__":
    main()
