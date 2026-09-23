#!/usr/bin/env python3
"""Split a genome FASTA into size-capped chunks for MINIPROT_ALIGN_CHUNK
(plan D: chunked miniprot) -- so each chunk's miniprot index is a fraction
of the whole genome's, capping peak alignment RAM (~10 GB per Gb of genome,
dominated by the index).

Greedily packs WHOLE sequences, in FASTA order, into chunk_000.fa,
chunk_001.fa, ...: a sequence joins the current chunk if doing so keeps that
chunk's cumulative length <= --chunk_bp, otherwise a new chunk starts. No
sequence is ever split across chunks, so a single sequence longer than
--chunk_bp still gets a whole chunk to itself (that chunk alone won't shrink
the index below one full sequence, but nothing else does either).

Two full passes over the input, both streaming: sequence data itself is
never held in memory (only the current line), which is what lets this scale
to a genome too large to buffer.
  1. record each sequence's id and length only (no sequence bytes kept), to
     decide the whole packing up front;
  2. re-read the FASTA and copy each sequence's header + lines straight into
     its already-decided chunk file.

Alongside the chunks: chunks.tsv (file<TAB>length_bp, one row per chunk, in
chunk order -- no header line, same convention as rename_sequences.py's
--out_sizes/--out_gaps) and total_length.txt (the whole genome's length, a
single integer -- MINIPROT_ALIGN_CHUNK needs this for its -G max-intron
formula, which is defined on the WHOLE genome's length, not a chunk's; see
that process and docs/plans/D_chunked_miniprot.md).
"""
import argparse
import gzip
import os
import sys

GZIP_MAGIC = b'\x1f\x8b'


def _opener_for(fasta_path):
    """Gzip is detected by magic bytes, not extension -- same convention as
    rename_sequences.py."""
    with open(fasta_path, 'rb') as f:
        magic = f.read(2)
    return gzip.open if magic == GZIP_MAGIC else open


def _header_id(line):
    """line is a stripped '>...' header (bytes); returns its first
    whitespace-delimited token, decoded, the same way rename_sequences.py
    and miniprot itself take a FASTA record's id."""
    header = line[1:].decode('utf-8', errors='replace')
    return header.split(None, 1)[0] if header else ''


def scan_lengths(fasta_path):
    """First pass: [(seq_id, length), ...] in FASTA order. Counts line
    lengths only -- no sequence bytes are kept."""
    opener = _opener_for(fasta_path)
    records = []
    seq_id = None
    length = 0
    with opener(fasta_path, 'rb') as f:
        for raw_line in f:
            line = raw_line.strip()
            if line[:1] == b'>':
                if seq_id is not None:
                    records.append((seq_id, length))
                seq_id = _header_id(line)
                length = 0
                continue
            if seq_id is None or not line:
                continue
            length += len(line)
        if seq_id is not None:
            records.append((seq_id, length))
    if not records:
        sys.exit(f"ERROR: no sequences found in {fasta_path}")
    return records


def assign_chunks(records, chunk_bp):
    """Greedy bin-packing of whole sequences, in FASTA order (see module
    docstring). Returns (seq_id -> chunk_index, [chunk_length_bp, ...])."""
    assignment = {}
    chunk_lengths = []
    current_chunk = -1
    current_len = 0
    for seq_id, length in records:
        if current_chunk == -1 or current_len + length > chunk_bp:
            current_chunk += 1
            chunk_lengths.append(0)
            current_len = 0
        assignment[seq_id] = current_chunk
        chunk_lengths[current_chunk] += length
        current_len += length
    return assignment, chunk_lengths


def write_chunks(fasta_path, assignment, n_chunks, outdir):
    """Second pass: stream each sequence's header + lines straight into its
    already-decided chunk file. Only the current line is ever held."""
    opener = _opener_for(fasta_path)
    paths = [os.path.join(outdir, f'chunk_{i:03d}.fa') for i in range(n_chunks)]
    outputs = [open(p, 'wb') for p in paths]
    try:
        current_out = None
        with opener(fasta_path, 'rb') as f:
            for raw_line in f:
                line = raw_line.strip()
                if line[:1] == b'>':
                    seq_id = _header_id(line)
                    current_out = outputs[assignment[seq_id]]
                    current_out.write(b'>')
                    current_out.write(line[1:])
                    current_out.write(b'\n')
                    continue
                if current_out is None or not line:
                    continue
                current_out.write(line)
                current_out.write(b'\n')
    finally:
        for out in outputs:
            out.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fasta', required=True, help='genome FASTA, optionally gzip-compressed')
    parser.add_argument('--chunk_bp', type=int, required=True,
                         help='target cumulative length (bp) per chunk; a single sequence '
                              'longer than this still gets a whole chunk to itself')
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    if args.chunk_bp <= 0:
        sys.exit(f"ERROR: --chunk_bp must be positive, got {args.chunk_bp}")

    os.makedirs(args.outdir, exist_ok=True)

    records = scan_lengths(args.fasta)
    assignment, chunk_lengths = assign_chunks(records, args.chunk_bp)
    n_chunks = len(chunk_lengths)
    write_chunks(args.fasta, assignment, n_chunks, args.outdir)

    total_length = sum(length for _seq_id, length in records)
    with open(os.path.join(args.outdir, 'chunks.tsv'), 'w') as f:
        for i, length in enumerate(chunk_lengths):
            f.write(f"chunk_{i:03d}.fa\t{length}\n")
    with open(os.path.join(args.outdir, 'total_length.txt'), 'w') as f:
        f.write(f"{total_length}\n")

    print(f"[split_genome] {len(records)} sequence(s), {total_length:,} bp total, "
          f"split into {n_chunks} chunk(s) (chunk_bp={args.chunk_bp:,})", file=sys.stderr)


if __name__ == '__main__':
    main()
