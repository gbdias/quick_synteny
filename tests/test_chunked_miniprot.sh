#!/usr/bin/env bash
# Acceptance checks 1-3 for plan D (chunked miniprot) --
# docs/plans/D_chunked_miniprot.md. Checks 4 (DAG -preview) and 5 (full
# local pipeline run, unchunked vs --miniprot_chunk_gb) are separate manual
# steps using the commands documented in that plan and in
# docs/plans/README.md's "Status" section.
#
# Needs docker (every miniprot call goes through the miniprot container, per
# docs/plans/README.md) and python3. Invoke as:
#   bash tests/test_chunked_miniprot.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MINIPROT_IMAGE='quay.io/biocontainers/miniprot:0.18--h577a1d6_0'
SEED=20260923
CHUNK_BP=500000   # forces the 3-chromosome mini_genome.fa fixture into >= 3 chunks

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/test_chunked_miniprot.XXXXXX")"
trap 'rm -rf "$WORKDIR"' EXIT
echo "[test] workdir: $WORKDIR"

n_checks=0
pass() { n_checks=$((n_checks + 1)); echo "[test] PASS: $1"; }
fail() { echo "[test] FAIL: $1" >&2; exit 1; }

miniprot() {
    docker run --rm --platform linux/amd64 \
        -v "$WORKDIR":"$WORKDIR" -w "$WORKDIR" \
        "$MINIPROT_IMAGE" miniprot "$@"
}

python3 "$REPO_ROOT/tests/fixtures/make_mini_genome.py" --seed "$SEED" --outdir "$WORKDIR/fixtures"
GENOME="$WORKDIR/fixtures/mini_genome.fa"
PROTEOME="$WORKDIR/fixtures/mini_proteome.faa"

# ---------------------------------------------------------------------------
# Check 1: -I vs an explicit -G computed by the plan's formula, on the WHOLE
# fixture genome, must give an identical GFF (after dropping ##PAF lines) --
# this is what confirms MINIPROT_ALIGN_CHUNK's -G replicates what -I would
# have done on the whole genome.
# ---------------------------------------------------------------------------
python3 "$REPO_ROOT/bin/split_genome.py" --fasta "$GENOME" --chunk_bp 999999999999 \
    --outdir "$WORKDIR/whole" >/dev/null
TOTAL_LENGTH="$(cat "$WORKDIR/whole/total_length.txt")"
G="$(python3 -c "
import math
L = $TOTAL_LENGTH
print(min(max(math.ceil(3.6 * math.sqrt(L)), 10000), 300000))
")"
echo "[test] genome length=$TOTAL_LENGTH -> G=$G"

miniprot -t 2 -I -N 5 --outs=0.7 --gff "$GENOME" "$PROTEOME" \
    > "$WORKDIR/whole_I.gff" 2> "$WORKDIR/whole_I.log"
miniprot -t 2 -G "$G" -N 5 --outs=0.7 --gff "$GENOME" "$PROTEOME" \
    > "$WORKDIR/whole_G.gff" 2> "$WORKDIR/whole_G.log"
grep -v '^##PAF' "$WORKDIR/whole_I.gff" > "$WORKDIR/whole_I.nopaf.gff"
grep -v '^##PAF' "$WORKDIR/whole_G.gff" > "$WORKDIR/whole_G.nopaf.gff"

if diff -q "$WORKDIR/whole_I.nopaf.gff" "$WORKDIR/whole_G.nopaf.gff" > /dev/null; then
    pass "check 1: -I and -G $G (formula) give byte-identical GFF on the whole genome"
else
    diff "$WORKDIR/whole_I.nopaf.gff" "$WORKDIR/whole_G.nopaf.gff" | head -20 || true
    fail "check 1: -I and explicit -G $G differ -- see diff above"
fi

# ---------------------------------------------------------------------------
# Check 2: chunked (>= 3 chunks) + merge vs unchunked, same -G, must give
# identical sets of (protein, chrom, start, end, strand, rank, score).
# ---------------------------------------------------------------------------
python3 "$REPO_ROOT/bin/split_genome.py" --fasta "$GENOME" --chunk_bp "$CHUNK_BP" \
    --outdir "$WORKDIR/chunks"
N_CHUNKS="$(wc -l < "$WORKDIR/chunks/chunks.tsv" | tr -d ' ')"
echo "[test] chunk_bp=$CHUNK_BP -> $N_CHUNKS chunk(s)"
if [ "$N_CHUNKS" -lt 3 ]; then
    fail "check 2: fixture only produced $N_CHUNKS chunk(s) at chunk_bp=$CHUNK_BP, need >= 3"
fi

CHUNK_GFFS=()
for CHUNK_FA in "$WORKDIR"/chunks/chunk_*.fa; do
    CHUNK_NAME="$(basename "$CHUNK_FA" .fa)"
    miniprot -t 2 -G "$G" -N 5 --outs=0.7 --gff "chunks/${CHUNK_NAME}.fa" "fixtures/$(basename "$PROTEOME")" \
        > "$WORKDIR/chunks/${CHUNK_NAME}.raw.gff" 2> "$WORKDIR/chunks/${CHUNK_NAME}.log"
    CHUNK_GFFS+=("$WORKDIR/chunks/${CHUNK_NAME}.raw.gff")
done

python3 "$REPO_ROOT/bin/merge_miniprot_gff.py" --gff "${CHUNK_GFFS[@]}" \
    --out "$WORKDIR/merged.gff" --outs 0.7 --max_secondary 5

if python3 - "$WORKDIR/whole_G.gff" "$WORKDIR/merged.gff" <<'PYEOF'
import re
import sys

ATTR_RE = re.compile(r'(\w+)=([^;]+)')


def extract(path):
    rows = set()
    with open(path) as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[2] != 'mRNA':
                continue
            attrs = dict(ATTR_RE.findall(fields[8]))
            protein = attrs['Target'].split()[0]
            rows.add((protein, fields[0], fields[3], fields[4], fields[6], attrs['Rank'], fields[5]))
    return rows


whole = extract(sys.argv[1])
chunked = extract(sys.argv[2])
if whole != chunked:
    print(f"MISMATCH: {len(whole - chunked)} hit(s) only in whole, "
          f"{len(chunked - whole)} only in chunked", file=sys.stderr)
    for row in sorted(whole - chunked)[:10]:
        print("  only in whole:", row, file=sys.stderr)
    for row in sorted(chunked - whole)[:10]:
        print("  only in chunked:", row, file=sys.stderr)
    sys.exit(1)
print(f"[test] identical: {len(whole)} hit(s) match exactly between whole-genome and chunked+merged")
PYEOF
then
    pass "check 2: chunked ($N_CHUNKS chunks) + merge == unchunked (same -G)"
else
    fail "check 2: chunked/unchunked (protein, chrom, start, end, strand, rank, score) sets differ (see above)"
fi

# ---------------------------------------------------------------------------
# Check 3: merge unit checks on a hand-made pair of chunk GFFs, covering:
#   - a chunk-level secondary that drops out under the GLOBAL outs filter
#     (proteinX: chunk_000's rank-2 hit, score 71, was kept locally against
#     chunk_000's own best of 100, but chunk_001 has a rank-1 hit of 150 for
#     the same protein, so the global best is 150 and 71 < 0.7*150);
#   - a protein exceeding 1 + max_secondary across chunks (proteinY has 5
#     candidates across the two chunks, 4 of which pass the outs filter;
#     with --max_secondary 2 only the top 3 may survive);
#   - ID collisions across chunks (both chunks independently number their
#     first mRNA ID=MP000001);
#   - Rank= rewritten on survivors (both mRNA and child CDS lines).
# ---------------------------------------------------------------------------
mkdir -p "$WORKDIR/merge_unit"
cat > "$WORKDIR/merge_unit/chunk_000.gff" <<'GFFEOF'
##gff-version 3
##PAF	proteinX	100	0	100	+	chrA	1000	10	110	100	100	0
chrA	miniprot	mRNA	11	110	100	+	.	ID=MP000001;Rank=1;Identity=1.0000;Positive=1.0000;Target=proteinX 1 100
chrA	miniprot	CDS	11	110	100	+	0	Parent=MP000001;Rank=1;Identity=1.0000;Target=proteinX 1 100
##PAF	proteinX	100	0	100	+	chrA	1000	200	271	71	71	0
chrA	miniprot	mRNA	201	271	71	+	.	ID=MP000002;Rank=2;Identity=0.9000;Positive=0.9500;Target=proteinX 1 71
chrA	miniprot	CDS	201	271	71	+	0	Parent=MP000002;Rank=2;Identity=0.9000;Target=proteinX 1 71
##PAF	proteinY	200	0	200	+	chrA	1000	300	500	90	90	0
chrA	miniprot	mRNA	301	500	90	+	.	ID=MP000003;Rank=1;Identity=1.0000;Positive=1.0000;Target=proteinY 1 200
chrA	miniprot	CDS	301	500	90	+	0	Parent=MP000003;Rank=1;Identity=1.0000;Target=proteinY 1 200
##PAF	proteinY	200	0	200	+	chrA	1000	600	700	80	80	0
chrA	miniprot	mRNA	601	700	80	+	.	ID=MP000004;Rank=2;Identity=0.9000;Positive=0.9200;Target=proteinY 1 200
chrA	miniprot	CDS	601	700	80	+	0	Parent=MP000004;Rank=2;Identity=0.9000;Target=proteinY 1 200
GFFEOF
cat > "$WORKDIR/merge_unit/chunk_001.gff" <<'GFFEOF'
##gff-version 3
##PAF	proteinX	100	0	100	-	chrB	1000	10	160	150	150	0
chrB	miniprot	mRNA	11	160	150	-	.	ID=MP000001;Rank=1;Identity=1.0000;Positive=1.0000;Target=proteinX 1 100
chrB	miniprot	CDS	11	160	150	-	0	Parent=MP000001;Rank=1;Identity=1.0000;Target=proteinX 1 100
##PAF	proteinY	200	0	200	-	chrB	1000	200	290	85	85	0
chrB	miniprot	mRNA	201	290	85	-	.	ID=MP000002;Rank=3;Identity=0.9000;Positive=0.9100;Target=proteinY 1 200
chrB	miniprot	CDS	201	290	85	-	0	Parent=MP000002;Rank=3;Identity=0.9000;Target=proteinY 1 200
##PAF	proteinY	200	0	200	-	chrB	1000	400	470	60	60	0
chrB	miniprot	mRNA	401	470	60	-	.	ID=MP000003;Rank=4;Identity=0.8000;Positive=0.8500;Target=proteinY 1 200
chrB	miniprot	CDS	401	470	60	-	0	Parent=MP000003;Rank=4;Identity=0.8000;Target=proteinY 1 200
##PAF	proteinY	200	0	200	-	chrB	1000	500	560	82	82	0
chrB	miniprot	mRNA	501	560	82	-	.	ID=MP000004;Rank=5;Identity=0.8500;Positive=0.9000;Target=proteinY 1 200
chrB	miniprot	CDS	501	560	82	-	0	Parent=MP000004;Rank=5;Identity=0.8500;Target=proteinY 1 200
GFFEOF

python3 "$REPO_ROOT/bin/merge_miniprot_gff.py" \
    --gff "$WORKDIR/merge_unit/chunk_000.gff" "$WORKDIR/merge_unit/chunk_001.gff" \
    --out "$WORKDIR/merge_unit/merged.gff" --outs 0.7 --max_secondary 2

if python3 - "$WORKDIR/merge_unit/merged.gff" <<'PYEOF'
import re
import sys

ATTR_RE = re.compile(r'(\w+)=([^;]+)')
path = sys.argv[1]

mrnas = []
children_rank = {}
with open(path) as f:
    for line in f:
        if not line or line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        attrs = dict(ATTR_RE.findall(fields[8]))
        if fields[2] == 'mRNA':
            mrnas.append((attrs['ID'], attrs['Target'].split()[0], attrs['Rank'], fields[5]))
        else:
            children_rank.setdefault(attrs['Parent'], []).append(attrs['Rank'])

errors = []

# -- chunk-level secondary dropped under the global outs filter: proteinX's
# only survivor must be chunk_001's hit (score 150); chunk_000's score-100
# and score-71 hits (kept locally, valid there) must both be gone.
x_hits = [m for m in mrnas if m[1] == 'proteinX']
if len(x_hits) != 1 or x_hits[0][0] != 'c1_MP000001':
    errors.append(f"proteinX: expected exactly [c1_MP000001] to survive, got {x_hits}")

# -- exceeding 1 + max_secondary across chunks: proteinY has 4 candidates
# passing --outs 0.7 (90, 85, 82, 80) but --max_secondary 2 caps at 3, so
# the score-80 hit (c0_MP000004) must be dropped even though it clears outs.
y_hits = sorted([m for m in mrnas if m[1] == 'proteinY'], key=lambda m: int(m[2]))
y_ids = [m[0] for m in y_hits]
if y_ids != ['c0_MP000003', 'c1_MP000002', 'c1_MP000004']:
    errors.append(f"proteinY: expected [c0_MP000003, c1_MP000002, c1_MP000004] survivors in rank order, got {y_ids}")

# -- ID collisions across chunks: both inputs had a bare ID=MP000001;
# prefixing must keep them distinguishable (c0_ / c1_), not silently merge.
all_ids = [m[0] for m in mrnas]
if len(all_ids) != len(set(all_ids)):
    errors.append(f"duplicate IDs in merged output: {all_ids}")
if 'c1_MP000001' not in all_ids:
    errors.append("expected the chunk_001-prefixed ID c1_MP000001 in the output")

# -- Rank= rewritten on mRNA AND children: c1_MP000002's original Rank was
# 3 (chunk-local); its new rank must be 2, and its CDS child must match.
c1_mp2 = next(m for m in mrnas if m[0] == 'c1_MP000002')
if c1_mp2[2] != '2':
    errors.append(f"c1_MP000002: expected rewritten Rank=2, got Rank={c1_mp2[2]}")
if children_rank.get('c1_MP000002') != ['2']:
    errors.append(f"c1_MP000002's CDS child: expected Rank=2, got {children_rank.get('c1_MP000002')}")

if errors:
    for e in errors:
        print("MISMATCH:", e, file=sys.stderr)
    sys.exit(1)
print(f"[test] all 4 merge unit properties hold ({len(mrnas)} surviving record(s))")
PYEOF
then
    pass "check 3: merge unit checks (chunk-local secondary drop, >1+max_secondary truncation, ID collisions, Rank= rewrite)"
else
    fail "check 3: merge unit checks failed (see above)"
fi

echo "[test] ALL $n_checks CHECKS PASSED"
