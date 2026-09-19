#!/usr/bin/env python3
"""Rename FASTA sequence IDs to short, informative, plot-friendly labels.

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
renamed FASTA, since the current bash pipeline has no such mapping and
renamed plot labels are otherwise untraceable back to source accessions.
"""
import argparse
import re
import sys

CHROMOSOME_RE = re.compile(r'chromosome[,:]?\s+([\w.]+)', re.IGNORECASE)
LINKAGE_GROUP_RE = re.compile(r'linkage\s+group[,:]?\s+([\w.]+)', re.IGNORECASE)
LG_SHORT_RE = re.compile(r'\bLG0*(\w+)\b', re.IGNORECASE)
UNLOCALIZED_RE = re.compile(r'unlocalized|unplaced', re.IGNORECASE)
MITOCHONDRION_RE = re.compile(r'mitochondri', re.IGNORECASE)
PLASTID_RE = re.compile(r'chloroplast|plastid', re.IGNORECASE)
ACCESSION_RE = re.compile(r'\.\d+$')  # e.g. NC_000913.3, CM012345.1


def parse_fasta_headers(fasta_path):
    """Return [(seq_id, description, length), ...] in file order, without
    loading full sequences into memory."""
    records = []
    seq_id = None
    description = ''
    length = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                if seq_id is not None:
                    records.append((seq_id, description, length))
                header = line[1:].rstrip('\n')
                parts = header.split(None, 1)
                seq_id = parts[0]
                description = parts[1] if len(parts) > 1 else ''
                length = 0
            else:
                length += len(line.strip())
        if seq_id is not None:
            records.append((seq_id, description, length))
    return records


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


def write_renamed_fasta(fasta_path, id_map, out_path):
    with open(fasta_path) as fin, open(out_path, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                old_id = line[1:].split(None, 1)[0]
                fout.write(f">{id_map[old_id]}\n")
            else:
                fout.write(line)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--out_fasta', required=True)
    parser.add_argument('--out_lookup', required=True)
    args = parser.parse_args()

    records = parse_fasta_headers(args.fasta)
    if not records:
        sys.exit(f"ERROR: no sequences found in {args.fasta}")

    assigned = assign_names(records)
    id_map = {old: new for old, new, _length, _type in assigned}

    write_renamed_fasta(args.fasta, id_map, args.out_fasta)

    with open(args.out_lookup, 'w') as f:
        f.write("original_id\tnew_id\tlength\ttype\n")
        for old_id, new_id, length, kind in assigned:
            f.write(f"{old_id}\t{new_id}\t{length}\t{kind}\n")

    print(f"[rename_sequences] renamed {len(assigned)} sequence(s): "
          f"{sum(1 for *_r, k in assigned if k.startswith('chromosome'))} chromosome, "
          f"{sum(1 for *_r, k in assigned if k == 'organelle')} organelle, "
          f"{sum(1 for *_r, k in assigned if k == 'preserved')} preserved, "
          f"{sum(1 for *_r, k in assigned if k == 'scaffold')} scaffold(fallback)",
          file=sys.stderr)


if __name__ == '__main__':
    main()
