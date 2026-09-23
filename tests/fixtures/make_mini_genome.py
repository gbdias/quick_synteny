#!/usr/bin/env python3
"""Deterministic test-fixture generator for the genome-scan tests
(tests/test_prepare_genome.py) and the chunked-miniprot tests
(tests/test_chunked_miniprot.sh). Standard library only, so it runs unmodified on both the
host's Python 3.10 and the pipeline containers' Python 3.13.

Two independent fixture sets, both written to --outdir:

  scan.fa / scan_crlf.fa / scan.fa.gz
      ~10 sequences exercising rename_sequences.py's header classification
      (chromosome, unlocalized, mitochondrion, a short preserved token, and
      several versioned accessions of different lengths for the scafN
      fallback order) and its N-run gap scan (runs shorter/equal/longer than
      typical --min_gap values, at sequence start/end, spanning one and
      several line wraps, an all-N sequence, and a mixed-case Nn run), plus
      formatting variants: line widths 60 and 80, no trailing newline at
      EOF, a CRLF copy, and a gzip copy. scan_crlf.fa and scan.fa.gz are
      byte-for-byte re-encodings of scan.fa's own content (CRLF line
      endings; gzip-compressed), not independently generated, so a scan of
      any of the three must agree exactly with a scan of the others.

  mini_genome.fa / mini_proteome.faa
      3 chromosomes of random background DNA carrying ~60 genes (each a
      random protein back-translated with random synonymous codons, split
      into 2-6 exons joined by GT...AG introns, placed on a random strand),
      ~10 paralog loci (a second, ~10%-diverged copy of one of the 60
      proteins, so miniprot's -N 5 turns up a same-protein secondary hit at
      a different locus), a few N-gaps in intergenic sequence, and a
      proteome where 5 of the 60 proteins also appear a second time as an
      N-terminally truncated "isoform" under a different accession (so the
      same genomic locus can get hit by two different query proteins).
"""
import argparse
import gzip
import os
import random
import sys

# ---------------------------------------------------------------------------
# scan.fa / scan_crlf.fa / scan.fa.gz
# ---------------------------------------------------------------------------

# (id, description, length, line_width, [(run_start, run_len, mixed_case), ...])
# min_seq_size 1000 keeps exactly chr3 (1200) and CM012345.1 (1050) -- every
# other record is deliberately < 1000bp so that filter has something to drop.
SCAN_RECORDS = [
    dict(id='NC_010001.1',
         desc='Example organism chromosome 3, whole genome shotgun sequence',
         length=1200, width=60,
         runs=[(0, 10, False),      # run at the very start, == min_gap(10)
               (300, 1, False)]),   # single-base run, == min_gap(1)
    dict(id='NW_020000002.1',
         desc='Example organism chromosome 3, unlocalized genomic scaffold',
         length=900, width=60,
         runs=[(300, 20, True),     # mixed Nn run, mid-sequence
               (800, 100, False)]), # run at the very end, == min_gap(100)
    dict(id='NC_010003.1',
         desc='Example organism mitochondrion, complete genome',
         length=1100, width=60,
         runs=[(50, 150, False),    # spans several line wraps (width 60)
               (595, 15, False)]),  # spans exactly one line wrap
    dict(id='chrZcustom',
         desc='Manually curated scaffold, not part of the NCBI submission',
         length=700, width=80,
         runs=[(300, 5, False)]),   # short run, clear of any wrap
    dict(id='CM012345.1',
         desc='Example organism unplaced genomic scaffold, whole genome shotgun sequence',
         length=1050, width=60,
         runs=[(500, 99, False)]),  # just under min_gap(100)
    dict(id='CM012346.2',
         desc='Example organism unplaced genomic scaffold, whole genome shotgun sequence',
         length=950, width=60,
         runs=[(400, 101, False)]), # just over min_gap(100)
    dict(id='JAAXYZ010000001.1',
         desc='Example organism unplaced genomic scaffold, whole genome shotgun sequence',
         length=850, width=60,
         runs=[]),                  # no N's at all
    dict(id='NW_999999.1',
         desc='Example organism unplaced genomic scaffold, whole genome shotgun sequence',
         length=650, width=80,
         runs=[(0, 90, False)]),    # run at the start, spanning a wrap (width 80)
    dict(id='NC_999999.10',
         desc='Example organism unplaced genomic scaffold, whole genome shotgun sequence',
         length=550, width=60,
         runs=[(545, 5, False)]),   # run reaching exactly to the end
    dict(id='CM099999.1',
         desc='Example organism unplaced genomic scaffold, whole genome shotgun sequence',
         length=400, width=60,
         runs=[(0, 400, False)]),   # entirely N
]


def _build_scan_content(rng):
    out = bytearray()
    for rec in SCAN_RECORDS:
        seq = bytearray(rng.choices(b'ACGT', k=rec['length']))
        for start, run_len, mixed in rec['runs']:
            if mixed:
                patt = bytes(78 if i % 2 == 0 else 110 for i in range(run_len))  # 'N'/'n'
            else:
                patt = b'N' * run_len
            seq[start:start + run_len] = patt
        out += f">{rec['id']} {rec['desc']}\n".encode('ascii')
        width = rec['width']
        for i in range(0, rec['length'], width):
            out += bytes(seq[i:i + width])
            out += b'\n'
    if out.endswith(b'\n'):
        del out[-1]  # last sequence's last line ends the file with no trailing newline
    return bytes(out)


def write_scan_fixtures(seed, outdir):
    rng = random.Random(seed)
    content = _build_scan_content(rng)
    with open(os.path.join(outdir, 'scan.fa'), 'wb') as f:
        f.write(content)
    with open(os.path.join(outdir, 'scan_crlf.fa'), 'wb') as f:
        f.write(content.replace(b'\n', b'\r\n'))
    with gzip.open(os.path.join(outdir, 'scan.fa.gz'), 'wb') as f:
        f.write(content)


# ---------------------------------------------------------------------------
# mini_genome.fa / mini_proteome.faa
# ---------------------------------------------------------------------------

CODON_TABLE = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAT', 'AAC'],
    'D': ['GAT', 'GAC'],
    'C': ['TGT', 'TGC'],
    'Q': ['CAA', 'CAG'],
    'E': ['GAA', 'GAG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'K': ['AAA', 'AAG'],
    'M': ['ATG'],
    'F': ['TTT', 'TTC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
}
AA_ALPHABET = list(CODON_TABLE.keys())

_COMPLEMENT = bytes.maketrans(b'ACGTNacgtn', b'TGCANtgcan')


def _revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def _random_protein(rng, length):
    # a leading Met is cosmetic (miniprot needs no start codon), just realistic
    return 'M' + ''.join(rng.choice(AA_ALPHABET) for _ in range(length - 1))


def _mutate_protein(rng, protein, rate):
    chars = list(protein)
    for i, aa in enumerate(chars):
        if rng.random() < rate:
            chars[i] = rng.choice([a for a in AA_ALPHABET if a != aa])
    return ''.join(chars)


def _back_translate(rng, protein):
    return ''.join(rng.choice(CODON_TABLE[aa]) for aa in protein).encode('ascii')


def _split_lengths(rng, total, n, min_len=15):
    """n positive part-lengths summing to total, each >= min_len when that's
    actually achievable -- exact exon boundaries don't matter, only that
    there are 2-6 of them and none is a degenerate sliver."""
    if n == 1:
        return [total]
    min_len = min(min_len, total // n)
    weights = [rng.uniform(0.6, 1.4) for _ in range(n)]
    total_weight = sum(weights)
    parts = [max(min_len, round(total * w / total_weight)) for w in weights]
    parts[-1] += total - sum(parts)
    if parts[-1] < min_len:
        biggest = max(range(n - 1), key=lambda i: parts[i])
        shortfall = min_len - parts[-1]
        parts[biggest] -= shortfall
        parts[-1] += shortfall
    return parts


def _make_intron(rng):
    length = rng.randint(80, 3000)
    middle = ''.join(rng.choice('ACGT') for _ in range(length - 4)).encode('ascii')
    return b'GT' + middle + b'AG'


def _build_gene_block(rng, protein):
    """Exons+introns for one gene, in coding (5'->3') orientation -- the
    caller reverse-complements this before placing a '-' strand gene on the
    genome's + strand."""
    cds = _back_translate(rng, protein)
    n_exons = rng.randint(2, 6)
    exon_lens = _split_lengths(rng, len(cds), n_exons)
    exons, pos = [], 0
    for elen in exon_lens:
        exons.append(cds[pos:pos + elen])
        pos += elen
    block = bytearray()
    for i, exon in enumerate(exons):
        block += exon
        if i < len(exons) - 1:
            block += _make_intron(rng)
    return bytes(block)


def _random_dna(rng, length):
    return bytes(rng.choices(b'ACGT', k=length)) if length > 0 else b''


def _build_chromosome(rng, target_length, tasks, n_gaps):
    seq = bytearray()
    placed = []
    intergenic_spans = []
    for task in tasks:
        gap_start = len(seq)
        seq += _random_dna(rng, rng.randint(500, 5000))
        intergenic_spans.append((gap_start, len(seq)))
        block = task['block'] if task['strand'] == '+' else _revcomp(task['block'])
        start = len(seq)
        seq += block
        placed.append({**task, 'start': start, 'end': len(seq)})
    if len(seq) < target_length:
        pad_start = len(seq)
        seq += _random_dna(rng, target_length - len(seq))
        intergenic_spans.append((pad_start, len(seq)))

    candidates = [s for s in intergenic_spans if s[1] - s[0] >= 300]
    rng.shuffle(candidates)
    for span_start, span_end in candidates[:n_gaps]:
        span_len = span_end - span_start
        n_len = min(rng.randint(100, 2000), span_len - 50)
        if n_len < 50:
            continue
        n_start = span_start + rng.randint(0, span_len - n_len)
        seq[n_start:n_start + n_len] = b'N' * n_len

    return bytes(seq), placed


def _write_wrapped(f, seq, width):
    for i in range(0, len(seq), width):
        f.write(seq[i:i + width])
        f.write(b'\n')


def generate_miniprot_fixture(seed, outdir, n_canonical=60, n_paralogs=10, n_isoforms=5):
    rng = random.Random(seed)

    canonical = []
    for i in range(1, n_canonical + 1):
        protein = _random_protein(rng, rng.randint(200, 500))
        canonical.append({'id': f'gene{i:04d}', 'protein': protein})

    tasks = []
    for gene in canonical:
        tasks.append({'id': gene['id'], 'kind': 'canonical',
                       'block': _build_gene_block(rng, gene['protein']),
                       'strand': rng.choice('+-')})
    for gene in rng.sample(canonical, n_paralogs):
        mutated = _mutate_protein(rng, gene['protein'], rate=0.10)
        tasks.append({'id': gene['id'], 'kind': 'paralog',
                       'block': _build_gene_block(rng, mutated),
                       'strand': rng.choice('+-')})
    rng.shuffle(tasks)

    n_chroms = 3
    chrom_tasks = [[] for _ in range(n_chroms)]
    for i, task in enumerate(tasks):
        chrom_tasks[i % n_chroms].append(task)

    placements = []
    with open(os.path.join(outdir, 'mini_genome.fa'), 'wb') as f:
        for ci in range(n_chroms):
            chrom_name = f'chr{ci + 1}'
            target_length = rng.randint(300_000, 800_000)
            seq, placed = _build_chromosome(rng, target_length, chrom_tasks[ci], n_gaps=rng.randint(1, 2))
            f.write(f'>{chrom_name}\n'.encode('ascii'))
            _write_wrapped(f, seq, width=70)
            for p in placed:
                placements.append({**p, 'chrom': chrom_name})

    isoform_genes = rng.sample(canonical, n_isoforms)
    with open(os.path.join(outdir, 'mini_proteome.faa'), 'w') as f:
        for gene in canonical:
            f.write(f">{gene['id']}\n")
            for i in range(0, len(gene['protein']), 60):
                f.write(gene['protein'][i:i + 60] + '\n')
        for gene in isoform_genes:
            protein = gene['protein']
            trunc_start = rng.randint(int(0.1 * len(protein)), int(0.3 * len(protein)))
            truncated = protein[trunc_start:]
            f.write(f">{gene['id']}_isoform2\n")
            for i in range(0, len(truncated), 60):
                f.write(truncated[i:i + 60] + '\n')

    return {
        'n_canonical': n_canonical,
        'n_paralogs': n_paralogs,
        'n_isoforms': n_isoforms,
        'canonical_ids': [g['id'] for g in canonical],
        'placements': placements,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    write_scan_fixtures(args.seed, args.outdir)
    info = generate_miniprot_fixture(args.seed + 1, args.outdir)
    print(f"[make_mini_genome] wrote scan.fa(+crlf/gz) and mini_genome.fa/mini_proteome.faa "
          f"({info['n_canonical']} genes, {info['n_paralogs']} paralogs, {info['n_isoforms']} isoforms) "
          f"to {args.outdir}", file=sys.stderr)


if __name__ == '__main__':
    main()
