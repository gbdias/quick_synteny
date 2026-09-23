#!/usr/bin/env python3
"""Checks for bin/extract_hits.py against a hand-made miniprot-style GFF
(docs/specs/hit_table.md section 1). Plain asserts; run with
`python3 tests/test_extract_hits.py`."""
import gzip
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, 'bin', 'extract_hits.py')


def mrna(chrom, start, end, strand, mid, target, rank=1, positive=0.95, identity=0.9, score=1000):
    return (f"{chrom}\tminiprot\tmRNA\t{start}\t{end}\t{score}\t{strand}\t.\t"
            f"ID={mid};Rank={rank};Identity={identity:.4f};Positive={positive:.4f};Target={target} 1 300\n")


def cds(chrom, start, end, strand, mid):
    return f"{chrom}\tminiprot\tCDS\t{start}\t{end}\t500\t{strand}\t0\tParent={mid};Rank=1;Identity=0.9\n"


GFF = ''.join([
    '##gff-version 3\n',
    '##PAF\tignored\n',
    # host gene on '+', two isoforms sharing exon 2 -> one locus
    mrna('chr1', 1001, 50000, '+', 'MP1', 'hostA'),
    cds('chr1', 1001, 1200, '+', 'MP1'), cds('chr1', 20001, 20300, '+', 'MP1'), cds('chr1', 49801, 50000, '+', 'MP1'),
    mrna('chr1', 20101, 50000, '+', 'MP2', 'hostB', positive=0.91),
    cds('chr1', 20101, 20300, '+', 'MP2'), cds('chr1', 49801, 50000, '+', 'MP2'),
    # gene nested in the host's first intron, opposite strand -> its own locus
    mrna('chr1', 5001, 9000, '-', 'MP3', 'nested'),
    cds('chr1', 5001, 5300, '-', 'MP3'), cds('chr1', 8701, 9000, '-', 'MP3'),
    # same-strand gene nested in the intron, no exon overlap -> its own locus too
    mrna('chr1', 10001, 12000, '+', 'MP4', 'nestedSame'),
    cds('chr1', 10001, 12000, '+', 'MP4'),
    # hit on a sequence dropped by --min_seq_size (absent from chrom.sizes)
    mrna('scaf99', 1, 900, '+', 'MP5', 'lost'),
    cds('scaf99', 1, 900, '+', 'MP5'),
    # second chromosome listed first in chrom.sizes; secondary hit; no CDS rows
    mrna('chr2', 301, 900, '-', 'MP6', 'hostA', rank=2, positive=0.61, score=700),
])

SIZES = 'chr2\t100000\nchr1\t200000\n'


def run(gff_text, sizes_text):
    with tempfile.TemporaryDirectory() as d:
        g, s, o = (os.path.join(d, n) for n in ('in.gff', 'chrom.sizes', 'out.tsv.gz'))
        with open(g, 'w') as f:
            f.write(gff_text)
        with open(s, 'w') as f:
            f.write(sizes_text)
        subprocess.run([sys.executable, SCRIPT, '--gff', g, '--chrom_sizes', s, '--out', o],
                       check=True, capture_output=True)
        with open(o, 'rb') as f:
            raw = f.read()
    return raw, gzip.decompress(raw).decode().splitlines()


def main():
    raw, lines = run(GFF, SIZES)
    assert lines[0] == 'chrom\tstart\tend\tstrand\tpositive\tidentity\tscore\trank\tprotein\tlocus', lines[0]
    rows = [l.split('\t') for l in lines[1:]]
    by_prot = {r[8]: r for r in rows if not (r[8] == 'hostA' and r[0] == 'chr2')}

    assert len(rows) == 5, rows                                   # scaf99 dropped
    assert all(r[0] != 'scaf99' for r in rows)
    assert [r[0] for r in rows] == ['chr2', 'chr1', 'chr1', 'chr1', 'chr1'], 'chrom.sizes order, then start'
    assert [int(r[1]) for r in rows[1:]] == [1000, 5000, 10000, 20100], 'GFF start - 1'
    assert by_prot['hostA'][2] == '50000' and by_prot['hostA'][6] == '1000'

    loc = {p: int(r[9]) for p, r in by_prot.items()}
    assert loc['hostA'] == loc['hostB'], 'isoforms sharing a CDS are one locus'
    assert len({loc['hostA'], loc['nested'], loc['nestedSame']}) == 3, 'nested genes stay separate'
    chr2 = [r for r in rows if r[0] == 'chr2'][0]
    assert chr2[9] == '0' and chr2[7] == '2' and chr2[4] == '0.6100', 'chr2 comes first in chrom.sizes'
    assert sorted(int(r[9]) for r in rows) == [0, 1, 1, 2, 3], 'loci numbered in genomic order'

    raw2, _ = run(GFF, SIZES)
    assert raw == raw2, 'byte-identical output for identical input'
    print('test_extract_hits: all checks passed')


if __name__ == '__main__':
    main()
