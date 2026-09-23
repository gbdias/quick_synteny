# quick_synteny — NAR Genomics and Bioinformatics Application Note plan

Target: **NAR Genomics and Bioinformatics**, Application Note.
Limits: max 3000 words, **3 figures/tables**, 30 references.
Emphasis: FAIR/open-source, reproducibility — plays to the existing
Nextflow + containers + pinned-version setup directly.

This supersedes `BENCHMARK_PROTOCOL.md` for the current submission target.
That document's full grid is shelved, not executed, kept only as a reference
for a possible future full benchmark paper.

---

## 1. The three figures/tables

Everything in the main text has to earn one of three slots. Current
assignment (superseding an earlier draft that put the memory curve in slot 1
as a figure and reserved a slot each for feature/quantitative comparison
tables):

| # | Content | Type | Status |
|---|---|---|---|
| 1 | Algorithm explanation (miniprot double-duty as ortholog finder, gap-bounded LMS chaining) | Figure | not started |
| 2 | Multipanel anatomy of the interactive result (ring / zoom / dotplot) | Figure | not started |
| 3 | Memory-vs-genome-size scaling, rice → axolotl, with "fits under 128 GB?" | Table | **done** — see §2 |

**Consequence:** with all three slots spent this way, there is no dedicated
slot left for a feature-comparison table (quick_synteny vs hobrac vs MCScanX
vs i-ADHoRe vs JCVI) or a quantitative hobrac/MCScanX benchmark table. Those
numbers go into prose instead — cited inline in Results/Discussion, with full
data available in the repo — rather than being cut outright. Flagged when
this restructure was requested; revisit if reviewers ask for a dedicated
comparison table and a slot needs to be freed up.

A homeolog ring-plot screenshot is strong material but doesn't need its own
slot — it can run as a sub-panel of Figure 1 or inline in the text alongside
Table 3's *suecica* row, keeping the count at exactly three.

---

## 2. Table — memory scaling

**What it is:** peak RSS vs. genome size across rice, maize, human,
*Ae. tauschii*, wheat, *Pleurodeles*, axolotl, plus a "fits under 128 GB?"
column against a common-HPC-node reference. This became the Table (not
Figure 1 — see §1 for the current figure assignments).

**Status: done.** Clean run, 2026-09-20, `node-c01` (384 GB), exclusive
allocation, all 7 genomes `OK`:

| genome | size (Gb) | peak RSS (GB) | GB/Gb | under 128 GB? |
|---|---|---|---|---|
| rice | 0.37 | 3.7 | 9.9 | yes |
| human | 3.30 | 26.8 | 8.1 | yes |
| maize | 2.18 | 25.3 | 11.6 | yes |
| *Ae. tauschii* | 4.22 | 50.4 | 11.9 | yes |
| wheat | 14.57 | 132.5 | 9.1 | **no** |
| *P. waltl* | 19.77 | 230.3 | 11.7 | **no** |
| axolotl | 28.21 | 260.6 | 9.2 | **no** |

Fit (least squares through the origin): **9.91 GB RAM per Gb genome**;
predicted 128 GB ceiling at **12.9 Gb**.

Two findings worth a sentence each in the text:
- `index_rss_gb` ≈ `align_rss_gb` for every genome (axolotl: 260.56 GB either
  way) — peak RAM is dominated by the genome index, not by aligning the
  proteome against it.
- Nothing actually OOM'd here. A **first attempt** on shared (non-exclusive)
  nodes produced nonsense across the board — rice OOM'd at 2.82 GB against a
  16 GB request — purely from node-level contention with co-scheduled jobs,
  not a real memory limit. That run was discarded outright, not corrected;
  see `BENCHMARK_PROTOCOL.md` prediction 8 for the full account. The earlier
  "ceiling around 15–20 Gb" estimate was an artifact of that contamination
  and is retracted, not revised.

**What's left:** the `-M` sweep (miniprot's k-mer sampling exponent, added as
a real pipeline flag — see `modules/local/miniprot_align.nf`), to see whether
it brings wheat/*Pleurodeles*/axolotl back under 128 GB (or further, toward
laptop-plausible RAM), and at what synteny-concordance cost. Tooling is built
(`benchmark/miniprot_m_sweep/`) but not yet run.

---

## 3. Table 1 — feature comparison

Purely qualitative, zero compute. Columns: annotation required (Y/N),
proteome/BUSCO anchor source, handles polyploidy/homeologs (Y/N), typical
runtime class, reference-genome selection (automatic/manual).

Rows: quick_synteny, hobrac, MCScanX, i-ADHoRe, JCVI/MCscan (the legacy
predecessor — also justifies the "replaces our own earlier pipeline" claim
with the ~35–40 min → ~2–3 min figure already in the README).

This table does the paper's structural work: it's what lets the text say
"quick_synteny is the only annotation-free tool that also handles
polyploidy" without needing a benchmark run to prove the *absence* of a
hobrac feature — that's a documented limitation in their own paper, citable
directly.

---

## 4. Table 2 — quantitative comparison

Reuses the two datasets already validated in the README (*A. thaliana* ×
*A. lyrata*; *A. thaliana* × *A. suecica*). No new dataset acquisition.

New runs needed, on data already staged:
- hobrac on both pairs (BUSCO lineage: brassicales)
- MCScanX on both pairs (annotation already exists for both — this is the
  ceiling row)

Columns: wall-clock, peak RSS, anchor count, block count, homeolog pairs
found (suecica row only — this is the cell MCScanX/hobrac cannot fill, and
that empty cell is doing real argumentative work).

This is roughly 4 tool runs total, not a grid — an afternoon of compute, not
a cluster campaign.

---

## 5. What gets cut from the shelved plan

Explicitly **not** happening for this submission (revisit only for a future
full benchmark paper, per the user's choice to skip the supplementary
stress-axis version too):

- The OFAT stress grid (divergence, contiguity, ploidy, BUSCO-representation,
  proteome-choice axes)
- Curated karyotype ground truth (Muller elements, YGOB, ALGs) and all
  simulation-based truth sets
- The 2×2 anchor/chaining ablation
- i-ADHoRe as a scored polyploid ceiling (mentioned qualitatively in Table 1
  only, not run)
- compleasm-backed hobrac re-benchmark
- Pre-registration

None of this is wasted: `BENCHMARK_PROTOCOL.md` keeps the full design for
whenever (if ever) a standalone benchmark paper is worth doing separately.

---

## 6. Manuscript skeleton (~3000 words)

1. **Introduction** (~400 words) — the gap: annotation-free, polyploid-aware,
   fast synteny. Cite hobrac as the direct peer, MCScanX/i-ADHoRe as the
   annotation-based standard, JCVI as the legacy approach this pipeline itself
   replaced.
2. **Implementation** (~1000 words) — miniprot double-duty as ortholog finder,
   gap-bounded LMS chaining, auto-tuned identity/anchor thresholds, taxonomy-
   guided comparison-genome discovery, homeolog detection mechanism, Nextflow/
   container architecture (the FAIR/reproducibility angle NARGAB weighs).
3. **Results** (~1200 words) — Table 1, Table 2, Figure 1, one paragraph per
   real-data validation, the ~35–40 min → ~2–3 min runtime story.
4. **Availability** (~100 words) — GitHub repo, license, container images,
   Nextflow version pin.
5. **Discussion / limitations** (~300 words) — the honest ones already in the
   README: divergent-pair false positive/negative rate vs. jcvi/MCScanX,
   miniprot's memory ceiling for very large genomes (Figure 1's own number),
   the ancestral-WGD confound in homeolog detection.

---

## 7. Immediate next steps

1. Run the three `sacct`/`.time`-file diagnostics on the existing memory-
   profile jobs (already given, not yet executed as far as this plan knows).
2. Re-submit only what's needed once the bug-fixed script's results are
   trusted.
3. Run hobrac + MCScanX once each on both existing datasets (Table 2).
4. Draft Table 1 (no compute — can start immediately, in parallel with 1–3).
