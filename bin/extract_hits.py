#!/usr/bin/env python3
"""Turn one genome's miniprot GFF into the per-genome hit table that
bin/chain.js chains -- the contract is docs/specs/hit_table.md section 1.

One row per miniprot mRNA (one alignment of one proteome protein) on a
sequence kept in chrom.sizes, plus a threshold-independent `locus` id:
hits on the same chrom and strand whose coding exons (CDS rows) overlap are
the same locus. Isoforms share coding exons, so they collapse into one
locus and can't count twice as independent synteny evidence; a gene nested
in another gene's intron shares no exon with it, so it stays separate --
which overlapping mRNA *spans* would not give, since in large genomes a
single span can cover several neighbouring genes.

Output is a gzip'd TSV written with a fixed gzip mtime, so identical input
gives byte-identical output.
"""
import argparse
import gzip
import re
import sys
from collections import defaultdict

ATTR_RE = re.compile(r'(\w+)=([^;]+)')
PARENT_RE = re.compile(r'(?:^|;)Parent=([^;]+)')
HEADER = ('chrom', 'start', 'end', 'strand', 'positive', 'identity', 'score', 'rank', 'protein', 'locus')


def read_chrom_sizes(path):
    names = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                names.append(line.split('\t')[0])
    return {name: i for i, name in enumerate(names)}


def parse_gff(path, chrom_index):
    """-> (hits, cds, n_dropped): hits are [chrom_idx, start, end, strand,
    positive, identity, score, rank, protein]; cds maps hit index -> list of
    (start, end) coding intervals."""
    hits, cds = [], defaultdict(list)
    id_to_hit = {}
    n_dropped = 0
    with open(path) as f:
        for line in f:
            if not line or line[0] == '#':
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            kind = fields[2]
            if kind == 'mRNA':
                ci = chrom_index.get(fields[0])
                if ci is None:
                    n_dropped += 1
                    continue
                attrs = dict(ATTR_RE.findall(fields[8]))
                protein = attrs.get('Target', '').split()
                if not protein:
                    continue
                identity = float(attrs['Identity']) if 'Identity' in attrs else None
                positive = float(attrs['Positive']) if 'Positive' in attrs else identity
                if positive is None:
                    positive = 0.0
                if identity is None:
                    identity = positive
                score = int(float(fields[5])) if fields[5] != '.' else 0
                id_to_hit[attrs.get('ID')] = len(hits)
                hits.append([ci, int(fields[3]) - 1, int(fields[4]), fields[6], positive, identity,
                             score, int(attrs.get('Rank', 1)), protein[0]])
            elif kind == 'CDS':
                parent = PARENT_RE.search(fields[8])
                h = id_to_hit.get(parent.group(1)) if parent else None
                if h is not None:
                    cds[h].append((int(fields[3]) - 1, int(fields[4])))
    return hits, cds, n_dropped


def assign_loci(hits, cds):
    """Connected components of same-(chrom, strand) hits with overlapping
    CDS intervals (an mRNA without CDS rows falls back to its own span).
    hits must already be in row order; returns one locus id per row,
    numbered by (chrom, locus min start, locus max end, smallest row)."""
    n = len(hits)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    by_strand = defaultdict(list)
    for i, h in enumerate(hits):
        for s, e in (cds.get(i) or [(h[1], h[2])]):
            by_strand[(h[0], h[3])].append((s, e, i))

    for intervals in by_strand.values():
        intervals.sort()
        head, reach = -1, -1
        for s, e, i in intervals:
            if s >= reach:
                head, reach = i, e
                continue
            a, b = find(head), find(i)
            if a != b:
                if b < a:
                    a, b = b, a
                parent[b] = a
            reach = max(reach, e)

    members = defaultdict(list)
    for i in range(n):
        members[find(i)].append(i)
    order = sorted(members.values(),
                   key=lambda rows: (hits[rows[0]][0], min(hits[r][1] for r in rows),
                                     max(hits[r][2] for r in rows), min(rows)))
    locus = [0] * n
    for lid, rows in enumerate(order):
        for r in rows:
            locus[r] = lid
    return locus, len(order)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gff', required=True, help='miniprot --gff output, sequence IDs already renamed')
    parser.add_argument('--chrom_sizes', required=True,
                        help='this genome\'s chrom.sizes -- hits on sequences not listed here are dropped')
    parser.add_argument('--out', required=True, help='output .hits.tsv.gz')
    args = parser.parse_args()

    chrom_index = read_chrom_sizes(args.chrom_sizes)
    chrom_names = sorted(chrom_index, key=chrom_index.get)
    hits, cds, n_dropped = parse_gff(args.gff, chrom_index)

    # row order (spec 1): chrom, start, end, protein, rank -- then score and
    # original position so the order is total even for degenerate duplicates
    order = sorted(range(len(hits)), key=lambda i: (hits[i][0], hits[i][1], hits[i][2], hits[i][8],
                                                    hits[i][7], -hits[i][6], i))
    rows = [hits[i] for i in order]
    row_cds = {new: cds[old] for new, old in enumerate(order) if old in cds}
    locus, n_loci = assign_loci(rows, row_cds)

    with open(args.out, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as gz:
        gz.write(('\t'.join(HEADER) + '\n').encode())
        for h, lid in zip(rows, locus):
            gz.write((f"{chrom_names[h[0]]}\t{h[1]}\t{h[2]}\t{h[3]}\t{h[4]:.4f}\t{h[5]:.4f}\t"
                      f"{h[6]}\t{h[7]}\t{h[8]}\t{lid}\n").encode())

    print(f"[extract_hits] {len(rows)} hit(s) on {len({h[0] for h in rows})} sequence(s) -> "
          f"{n_loci} locus/loci; {n_dropped} hit(s) on sequences absent from chrom.sizes dropped "
          f"-> {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
