#!/usr/bin/env python3
"""Rename FASTA sequence IDs to short, informative, plot-friendly labels, and
compute the two things everything downstream needs about the genome (its
sequence sizes and its assembly-gap runs) in the same single pass.

NCBI headers vary a lot: some carry an explicit chromosome/linkage-group
label in the description, some are bare accessions with nothing informative.
Priority order:
  1. explicit "chromosome N" / "linkage group N" / "LGN" in the header ->
     chrN (unlocalized/unplaced scaffolds tied to a chromosome get a
     _u<n> suffix, matching the original awk script's convention).
  2. mitochondrion -> chrM, chloroplast/plastid -> chrPltd.
  3. an ID that already looks like a short, clean token (not a versioned
     GenBank/RefSeq accession) is preserved as-is.
  4. everything else is sorted by sequence length descending and assigned
     scafN in that order.

Always writes a lookup TSV (original_id, new_id, length, type) alongside the
size/gap tables, since the current pipeline has no other mapping and renamed
plot labels are otherwise untraceable back to source accessions.

No FASTA is rewritten. Earlier versions of this pipeline wrote a renamed copy
of the genome, then read that copy twice more (faSize, then a separate
assembly-gap scan) -- three full passes over data that can be tens of GB.
This script reads the ORIGINAL fasta exactly once and emits the lookup plus
the size table (--out_sizes, replacing faSize -detailed | awk) and the gap
table (--out_gaps, replacing the old find_assembly_gaps.py) directly, keyed
by the renamed IDs; nothing downstream needs a renamed FASTA at all -- the
GFF from aligning against the ORIGINAL fasta is renamed in place instead
(see rename_gff.py).
"""
import argparse
import gzip
import re
import sys

CHROMOSOME_RE = re.compile(r'chromosome[,:]?\s+([\w.]+)', re.IGNORECASE)
LINKAGE_GROUP_RE = re.compile(r'linkage\s+group[,:]?\s+([\w.]+)', re.IGNORECASE)
LG_SHORT_RE = re.compile(r'\bLG0*(\w+)\b', re.IGNORECASE)
UNLOCALIZED_RE = re.compile(r'unlocalized|unplaced', re.IGNORECASE)
MITOCHONDRION_RE = re.compile(r'mitochondri', re.IGNORECASE)
PLASTID_RE = re.compile(r'chloroplast|plastid', re.IGNORECASE)
ACCESSION_RE = re.compile(r'\.\d+$')  # e.g. NC_000913.3, CM012345.1

GZIP_MAGIC = b'\x1f\x8b'
N_RUN_RE = re.compile(rb'[Nn]+')


def _opener_for(fasta_path):
    """Gzip is detected by magic bytes, not extension -- a caller may hand us
    a gzipped file under any name."""
    with open(fasta_path, 'rb') as f:
        magic = f.read(2)
    return gzip.open if magic == GZIP_MAGIC else open


def scan_fasta(fasta_path):
    """The only full read of the genome. Returns:
      records: [(seq_id, description, length), ...] in file order
      gaps: {seq_id: [(start, end), ...]} in position order, using the
            ORIGINAL sequence ids (renamed by the caller after assign_names())

    Gap semantics match the old find_assembly_gaps.py exactly: 0-based,
    half-open runs of [Nn] of any length (the --min_gap floor is applied by
    the caller at output time), contiguous across line wraps, measured on
    each line stripped of surrounding whitespace -- same as the old script's
    `line.strip()`, which is also what makes a CRLF input transparent here
    without any special-casing (bytes.strip() drops a trailing b'\\r' same
    as text-mode universal-newline translation would have already dropped
    it for the old script).

    Reads in binary mode (required to detect gzip by magic bytes rather than
    extension, and to scan for N-runs without a text-decoding detour); only
    header lines are decoded, as UTF-8 with errors='replace', since only
    those are ever split into an id/description.
    """
    opener = _opener_for(fasta_path)
    records = []
    gaps = {}

    seq_id = None
    description = ''
    length = 0
    run_start = None  # start offset (within the CURRENT sequence) of an open N-run, or None

    def close_current():
        if run_start is not None:
            gaps[seq_id].append((run_start, length))
        if seq_id is not None:
            records.append((seq_id, description, length))

    with opener(fasta_path, 'rb') as f:
        for raw_line in f:
            line = raw_line.strip()
            if line[:1] == b'>':
                close_current()
                header = line[1:].decode('utf-8', errors='replace')
                parts = header.split(None, 1)
                seq_id = parts[0] if parts else ''
                description = parts[1] if len(parts) > 1 else ''
                length = 0
                run_start = None
                gaps[seq_id] = []
                continue
            if seq_id is None or not line:
                continue

            line_len = len(line)
            # Most lines have no N at all -- skip the regex scan entirely for
            # those (this `in` check is what gets this past 100 MB/s; a
            # per-character Python loop over a whole chromosome is what took
            # the old find_assembly_gaps.py 57 minutes on axolotl).
            if b'N' in line or b'n' in line:
                matches = list(N_RUN_RE.finditer(line))
                if run_start is not None and (not matches or matches[0].start() != 0):
                    # the run open from a previous line didn't reach this one
                    gaps[seq_id].append((run_start, length))
                    run_start = None
                for i, m in enumerate(matches):
                    s, e = m.span()
                    if i == 0 and run_start is not None and s == 0:
                        if e == line_len:
                            pass  # still open, carries into the next line
                        else:
                            gaps[seq_id].append((run_start, length + e))
                            run_start = None
                        continue
                    if e == line_len:
                        run_start = length + s  # open, carries into the next line
                    else:
                        gaps[seq_id].append((length + s, length + e))
            elif run_start is not None:
                gaps[seq_id].append((run_start, length))
                run_start = None

            length += line_len
        close_current()

    return records, gaps


def classify(seq_id, description):
    """Return ('chromosome'|'organelle'|'preserved'|None, label) --
    None means "no confident label, fall back to size-based scafN"."""
    text = f"{seq_id} {description}"

    m = CHROMOSOME_RE.search(text) or LINKAGE_GROUP_RE.search(text) or LG_SHORT_RE.search(text)
    if m:
        token = m.group(1).rstrip('.,;')
        label = f"chr{token}"
        if UNLOCALIZED_RE.search(text):
            return ('chromosome_unplaced', label)
        return ('chromosome', label)

    if MITOCHONDRION_RE.search(text):
        return ('organelle', 'chrM')
    if PLASTID_RE.search(text):
        return ('organelle', 'chrPltd')

    if len(seq_id) <= 15 and not ACCESSION_RE.search(seq_id):
        return ('preserved', seq_id)

    return (None, None)


def assign_names(records):
    """records: [(seq_id, description, length), ...] -> [(seq_id, new_id, length, type), ...]"""
    assigned = []
    unplaced_counts = {}
    used_names = set()
    fallback = []

    for seq_id, description, length in records:
        kind, label = classify(seq_id, description)
        if kind == 'chromosome_unplaced':
            unplaced_counts[label] = unplaced_counts.get(label, 0) + 1
            new_id = f"{label}_u{unplaced_counts[label]}"
            assigned.append((seq_id, new_id, length, 'chromosome_unplaced'))
            used_names.add(new_id)
        elif kind in ('chromosome', 'organelle', 'preserved'):
            new_id = label
            if new_id in used_names:
                new_id = f"{new_id}_{seq_id}"
            assigned.append((seq_id, new_id, length, kind))
            used_names.add(new_id)
        else:
            fallback.append((seq_id, description, length))

    fallback.sort(key=lambda r: r[2], reverse=True)
    for i, (seq_id, description, length) in enumerate(fallback, start=1):
        new_id = f"scaf{i}"
        assigned.append((seq_id, new_id, length, 'scaffold'))

    order = {seq_id: i for i, (seq_id, *_r) in enumerate(records)}
    assigned.sort(key=lambda r: order[r[0]])
    return assigned


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fasta', required=True, help='genome FASTA, optionally gzip-compressed')
    parser.add_argument('--out_lookup', required=True)
    parser.add_argument('--out_sizes', required=True,
                         help='renamed_id\\tlength, natural FASTA order, filtered by --min_seq_size')
    parser.add_argument('--min_seq_size', type=int, required=True)
    parser.add_argument('--out_gaps', required=True,
                         help='renamed_id\\tstart\\tend for every sequence (not filtered by --min_seq_size)')
    parser.add_argument('--min_gap', type=int, required=True,
                         help='minimum run length (bp) of consecutive Ns to count as an assembly gap')
    args = parser.parse_args()

    records, gaps = scan_fasta(args.fasta)
    if not records:
        sys.exit(f"ERROR: no sequences found in {args.fasta}")

    assigned = assign_names(records)
    id_map = {old: new for old, new, _length, _type in assigned}

    with open(args.out_lookup, 'w') as f:
        f.write("original_id\tnew_id\tlength\ttype\n")
        for old_id, new_id, length, kind in assigned:
            f.write(f"{old_id}\t{new_id}\t{length}\t{kind}\n")

    with open(args.out_sizes, 'w') as f:
        for seq_id, _description, length in records:
            if length >= args.min_seq_size:
                f.write(f"{id_map[seq_id]}\t{length}\n")

    gap_count = 0
    gap_bp = 0
    with open(args.out_gaps, 'w') as f:
        for seq_id, _description, _length in records:
            for start, end in gaps[seq_id]:
                if end - start >= args.min_gap:
                    f.write(f"{id_map[seq_id]}\t{start}\t{end}\n")
                    gap_count += 1
                    gap_bp += end - start

    print(f"[rename_sequences] renamed {len(assigned)} sequence(s): "
          f"{sum(1 for *_r, k in assigned if k.startswith('chromosome'))} chromosome, "
          f"{sum(1 for *_r, k in assigned if k == 'organelle')} organelle, "
          f"{sum(1 for *_r, k in assigned if k == 'preserved')} preserved, "
          f"{sum(1 for *_r, k in assigned if k == 'scaffold')} scaffold(fallback)",
          file=sys.stderr)
    print(f"[rename_sequences] {gap_count} gap(s) >= {args.min_gap} bp "
          f"({gap_bp:,} bp total) written to {args.out_gaps}", file=sys.stderr)


if __name__ == '__main__':
    main()
