#!/usr/bin/env python3
"""Find assembly gaps (runs of N/n) in a FASTA, one row per run of at least
--min_gap consecutive N's.

Runs against the RENAMED genome FASTA (chr1, chr2, ... -- see
rename_sequences.py), not the original download, so gap coordinates land on
the same sequence names as chrom.sizes/links.tsv/the plot itself use
everywhere else.

Output TSV (chrom, start, end -- no header, one row per gap): 0-based,
half-open (BED-like), matching this pipeline's own links.tsv convention
(see build_synteny_blocks.py's parse_gff: GFF start - 1, GFF end unchanged).
"""
import argparse
import sys


def find_gaps(fasta_path, min_gap):
    """Streams the FASTA one line at a time (never holds a whole sequence in
    memory) tracking the current run of N/n bases; a run is flushed as a gap
    the moment it's broken by a non-N base or the sequence ends, so a gap
    that happens to straddle a FASTA line wrap is still counted as one
    contiguous run, not split at the wrap."""
    chrom = None
    pos = 0
    run_start = None

    def flush(gaps):
        if run_start is not None and pos - run_start >= min_gap:
            gaps.append((chrom, run_start, pos))

    gaps = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                flush(gaps)
                chrom = line[1:].split(None, 1)[0]
                pos = 0
                run_start = None
                continue
            for ch in line.strip():
                is_n = ch in ('N', 'n')
                if is_n and run_start is None:
                    run_start = pos
                elif not is_n and run_start is not None:
                    if pos - run_start >= min_gap:
                        gaps.append((chrom, run_start, pos))
                    run_start = None
                pos += 1
        flush(gaps)
    return gaps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--min_gap', type=int, default=100,
                         help='minimum run length (bp) of consecutive Ns to count as an '
                              'assembly gap (default: 100)')
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    gaps = find_gaps(args.fasta, args.min_gap)

    with open(args.out, 'w') as f:
        for chrom, start, end in gaps:
            f.write(f"{chrom}\t{start}\t{end}\n")

    total_bp = sum(end - start for _, start, end in gaps)
    print(f"[find_assembly_gaps] {len(gaps)} gap(s) >= {args.min_gap} bp "
          f"({total_bp:,} bp total) written to {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
