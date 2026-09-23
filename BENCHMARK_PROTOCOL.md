# Protein-based synteny benchmark — protocol

Draft protocol for a benchmark of protein/gene-anchored synteny tools, with
`quick_synteny` as one entrant. Nucleotide-only tools (minimap2, MUMmer,
D-GENIES, SyRI) are deliberately out of scope: they fail for reasons unrelated
to the design choices under test.

Status: **shelved for now.** The immediate goal changed to a NAR Genomics and
Bioinformatics Application Note (max 3000 words, 3 figures/tables) rather than
a standalone benchmark paper — see `APPLICATION_NOTE_PLAN.md` for the
right-sized package actually being built. This document is kept as a
reference design for a possible future full benchmark paper; nothing below is
being executed against the current submission target.

---

## 0. The question

Not "does quick_synteny win" — it will not win on accuracy against tools that
are handed an annotation, and the README already concedes that. The question
the grid answers is:

> Among tools that take **only an unannotated assembly**, how well is synteny
> recovered across divergence, contiguity, ploidy and genome size — and how far
> is that from what an annotation-based tool achieves on the same data?

The scored comparison is **restricted to annotation-free tools**, on matched
inputs. Annotation-based tools appear as a **reference ceiling**, not as
competitors (§1). Predictions are recorded in §11 **before** the grid runs.

---

## 1. Entrants

### Scored entrants — annotation-free, matched inputs

| Entrant | Anchors | Chaining |
|---|---|---|
| `quick_synteny` | miniprot, full proteome, `-N 5 --outs=0.7` | gap-bounded LMS |
| `hobrac` | BUSCO complete single-copy orthologues | Fisher exact + chromosome chains |

That is the entire annotation-free field: these two tools plus DIY BUSCO
scripting. **Say so explicitly in the paper** — n=2 is a fact about the field,
not a gap in the benchmark.

**Framing consequence.** With only two scored entrants, this is not a tool
bake-off and should not be written as one. Frame it as *characterising two
approaches to annotation-free synteny across conditions* — the contribution is
the **axes** (§5), not the ranking. A two-tool comparison across six stress
axes against a reference ceiling is a substantial and defensible paper; a
two-tool comparison on one dataset is not.

### Reference ceiling — annotation-based, NOT scored against

Run on the same datasets **with** their annotation, plotted as a labelled
reference line: *"upper bound; different input budget; not a like-for-like
comparison."*

| Tool | Role | Axes it runs on |
|---|---|---|
| MCScanX | the accuracy ceiling | baseline, divergence, annotation-quality |
| i-ADHoRe | the polyploid ceiling | baseline, ploidy |
| JCVI / MCscan | predecessor — runtime story only | baseline only |

Why keep them, having restricted the comparison:

1. **Scale.** "quick_synteny 0.78, hobrac 0.71" is uninterpretable without a
   ceiling — a reviewer cannot tell whether both tools are good or both bad.
2. **The polyploid claim gets much stronger.** Without i-ADHoRe it is "we do
   homeologs, hobrac doesn't" — beating a tool at something it never
   attempted. With it: *"we approach the annotation-based polyploid specialist
   without requiring annotation."*
3. **Pre-empts "you avoided the standard."** Legitimate scoping still reads as
   evasion if the standard is simply absent.

Dropped: GENESPACE (heavy OrthoFinder dependency, redundant with MCScanX once
Tier B is only a ceiling), WGDI (plant-only).

### Explicitly excluded

Nucleotide-only aligners; annotation-transfer tools (Liftoff, GeMoMa) — not
synteny callers; visualization-only tools (Circos, NGenomeSyn, SynVisio).

### Datasets must still be annotated

Every dataset is one where annotation *exists*, so the ceiling can be computed
— with the annotation then **withheld** from the scored entrants.

---

## 2. Ground truth

Consensus-among-tools is circular and reviewers know it. Use curated karyotype
systems as primary truth, simulation for the axes real data cannot supply.

### 2a. Curated real-data truth

| System | Taxon | Why it qualifies |
|---|---|---|
| **Muller elements** | *Drosophila* | 6 elements, correspondence known with near-certainty across the genus; deep divergence |
| **YGOB** | *Saccharomyces* | manually curated synteny across the yeast WGD; rare and excellent |
| **Merian elements** | Lepidoptera | curated ALG system |
| **Metazoan ALGs** | bivalves / molluscs | hobrac's own currency — cannot be accused of friendly ground |
| **Brassica 3:1 triplication** | Brassicaceae | well characterised, and a polyploidy signal |
| **Ensembl Compara synteny** | human/mouse | deeply curated, at scale |

### 2b. Simulation truth

For breakpoint-level and homeolog-level truth no real dataset provides:

- **Rearrangement simulation** — one annotated genome, a known set of
  inversions / translocations / fusions / fissions applied, annotation lifted
  over. Exact breakpoint truth.
- **Polyploid simulation** — duplicate, diverge under a substitution model,
  optionally fractionate. Yields **ground-truth homeologs**, the only way to
  prove the claim hobrac cannot contest.
- Seeds recorded. Replicate count set empirically — see §7.5.

Simulated truth sets get a Zenodo DOI; they are not re-derivable without the
exact seed and code version.

---

## 3. Metrics

**Score at the gene-pair (anchor) level, not the block level.** Block
boundaries are definitional and differ per tool; block-level comparison
penalises tools for granularity rather than correctness.

| Metric | Level | Applies to |
|---|---|---|
| Precision / recall / F1 on syntenic gene pairs | anchor | all |
| Chromosome-pair confusion matrix | chromosome | all (coarse, boundary-free, robust) |
| Block count, block N50, genome coverage in blocks | block | descriptive only, not scored |
| Rearrangement recall + breakpoint localisation error (bp) | breakpoint | simulation only |
| Homeolog chromosome-pair precision / recall | chromosome | ploidy axis |
| Wall-clock, CPU-hours, peak RSS, peak disk | resource | all |

Resource metrics come from Nextflow's `trace` (already enabled in
`quick_synteny`). Nextflow's `peak_rss` is `/proc`-polled and can undercount
short spikes — additionally wrap each tool in `/usr/bin/time -v` and prefer
that figure for the paper.

---

## 4. The common anchor format

One format, one scoring implementation. Every tool gets a thin adapter; the
scorer never sees tool-native output. This is the main fairness guarantee.

```
dataset_id  tool  replicate  gene_id  q_chrom  q_start  q_end
            s_chrom  s_start  s_end  block_id  score  status
```

`status` carries `OK` / `OOM` / `TIMEOUT` / `FAILED` so failures flow into the
results table rather than crashing the grid (§8).

---

## 5. Grid design — one factor at a time, NOT full factorial

Define one **baseline cell** and vary one axis at a time from it.

**Baseline:** *A. thaliana* × *A. lyrata*, chromosome-level, diploid, full
annotation, default parameters.

| Axis | Levels | Scored | Ceiling | Runs |
|---|---|---|---|---|
| baseline | 1 | 2 | 3 tools | 5 |
| **divergence** | genus → family → order → class | 2 | MCScanX | 12 |
| **contiguity** | chr-level, 10 Mb, 1 Mb, 100 kb N50 | 2 | — | 8 |
| **ploidy** | diploid, allotetraploid, hexaploid | 2 | i-ADHoRe | 9 |
| **annotation quality** | 100 / 75 / 50 % genes retained | — | MCScanX | 3 |
| **BUSCO representation** | well- vs poorly-represented phylum | 2 | — | 4 |
| **proteome choice** | 5 proteomes at varying distance | quick_synteny | — | 5 |
| **genome size** | 3 rungs, §6 | 2 | — | 6 |
| **simulation** | ~6 conditions × n seeds | 2 | — | ~24 |

≈ **76 runs.**

The annotation-quality axis runs *only* the ceiling tool — its purpose is to
price the hidden cost of requiring annotation, a statement about Tier B, not
about the scored entrants.

Axis rationale — each targets a stated weakness of one entrant:
- *contiguity*: hobrac lists it as a limitation
- *BUSCO representation*: hobrac concedes resolution depends on marker
  availability
- *annotation quality*: prices Tier B's input requirement
- *proteome choice*: `quick_synteny`'s **own** biggest methodological
  vulnerability — test it before a reviewer does

Since the contribution is the axes rather than a ranking, these six are the
paper's backbone and should each get a figure.

---

## 6. Genome-size ladder

Three rungs, spanning two and a half orders of magnitude.

| # | Pair | Sizes | Projected miniprot peak RSS | On 128 GB nodes |
|---|---|---|---|---|
| 1 | *A. thaliana* × *A. lyrata* | 0.12 / 0.21 Gb | ~2 GB | trivial |
| 2 | human × mouse | 3.1 / 2.7 Gb | **~22 GB (measured)** | fine; curated truth at scale |
| 3 | *P. waltl* × *A. mexicanum* | 20.3 / 32 Gb | ~140 / ~220 GB | **expected to exceed 128 GB** |

Projection basis: miniprot measured at 21.8–22.5 GB peak for 25,007 proteins
against the 3.1 Gb human genome → ~7 GB per Gb of genome, index-dominated.
**Measure, do not assume** — miniprot may batch internally and beat the linear
extrapolation.

Notes:
- Rung 3 pairs axolotl with *Pleurodeles* because axolotl has no sensible
  partner alone; both are chromosome-scale urodeles. Divergence is deep
  (~150 Mya), additionally testing protein-anchor divergence tolerance at
  scale.
- **hobrac is expected to win rung 3**: BUSCO is slow but memory-bounded and
  chunked, so it will complete where miniprot OOMs. Record this honestly.

### The gap between rung 2 and rung 3

At three rungs the ladder jumps from 22 GB (comfortable) straight to 140+ GB
(fails), with nothing near the 128 GB boundary. That gives a pass/fail result
but no located ceiling — "it works at 3 Gb and fails at 20 Gb" is a weak
statement.

**Fix it cheaply with an index-only memory profile.** Add a tiny sub-workflow
that runs *only* `miniprot -d` (index construction, no alignment, no BUSCO, no
scoring) across a size series — e.g. rice 0.4 Gb, maize 2.3 Gb, human 3.1 Gb,
*Ae. tauschii* 4.3 Gb, wheat 16 Gb, *P. waltl* 20.3 Gb, axolotl 32 Gb.

Each is minutes of CPU and one process. The output is a proper RAM-vs-genome-
size curve that **locates the ceiling to within a Gb or two**, rather than
bracketing it between 3 and 20. This is the highest value-per-CPU-hour item in
the whole protocol.

### Mitigation if miniprot OOMs — and a real caveat

The obvious fix is chunking the genome and merging coordinates. **But**:
miniprot's `-N 5` secondary-hit ranking is computed per invocation, so
chunking makes it per-chunk. `quick_synteny`'s homeolog detection depends on
one protein hitting chr9 *and* chr10 **in the same alignment run** — if
homeologous chromosomes land in different chunks, the homeolog signal is
destroyed.

Options, in order of preference:
1. Chunk, retain all chunks' hits, **re-rank globally** before chaining.
2. Chunk by homeolog-group where known — fragile, requires prior knowledge.
3. Accept the ceiling, document it as a hard limitation.

Option 1 is probably a feature `quick_synteny` needs regardless.

---

## 7. Making it cheaper

### 7.1 Cache BUSCO and miniprot across the whole grid — biggest single win

BUSCO is ~90 % of hobrac's CPU time (7.91 of 8.46 CPU-hours in their own
published run). The **same genome recurs in 10–15 grid cells**, and BUSCO's
result depends only on `(genome, lineage dataset)` — not on which axis is
being varied.

Run it **once per (genome, lineage)** and reuse. Same for miniprot GFFs, keyed
on `(proteome, genome, params)`, and for MCScanX's all-vs-all BLASTP.

Mechanism: Nextflow `storeDir` on a shared cache directory, keyed by a content
hash of the inputs. `storeDir` — not `-resume` — because it survives across
separate pipeline invocations and across the different stress workflows.

hobrac accepts pre-computed BUSCO results, so this needs no tool modification.

### 7.2 Subset chromosomes on the stress axes

The divergence, contiguity, annotation-quality and BUSCO-representation axes
measure *relative* differences between tools. They do not need whole genomes —
a consistent subset of ~3 mutually syntenic chromosomes per genome gives a
valid relative comparison at 3–5× lower cost.

Reserve whole genomes for the size ladder and the headline baseline numbers.
Never compare a subset figure against a whole-genome figure; mark subset runs
in the tidy TSV.

### 7.3 Stage the grid in three phases

The failure mode to avoid is burning hundreds of CPU-hours and then finding
the normalizer mis-parsed one tool's output.

1. **Phase 1** — the baseline cell only. Cheap, and it exercises every module
   end-to-end including scoring and collation.
2. **Phase 2** — all stress axes, small genomes, subset chromosomes.
3. **Phase 3** — rungs 2–3 and the index-only memory profile. Expensive, run
   last, only once everything else is stable.

### 7.4 Decide rung 3 in two stages

**Rung 3 is now ~60–70 % of the entire compute budget** — BUSCO on a 20 Gb and
a 32 Gb genome is on the order of 500–1,000 CPU-hours, for what ends up as a
handful of table cells.

Split the decision:

1. Run `quick_synteny` first. If it OOMs, that costs ~1 CPU-hour (the index
   build dies early) and your own ceiling is established.
2. *Then* decide whether hobrac's full run is worth ~600 CPU-hours for the
   "they complete where we fail" cell — or whether citing their published
   scaling figures suffices.

Running it is more convincing; citing is defensible. Make the call with the
rung-1/2 numbers already in hand.

### 7.5 Set the replicate count empirically

Run a pilot with 5 seeds on 2 simulation conditions, measure between-seed
variance, and justify n=2 or n=3 from that. Most of these pipelines are
deterministic given a genome, so n=2 is likely defensible and halves the
block. Trim to ~6 conditions.

### 7.6 Cluster hygiene

- `process.scratch = true` — node-local disk, avoids NFS thrash on 20–32 Gb
  genomes. Matters enormously at rung 3.
- Freeze datasets once into a shared read-only directory, checksummed. Never
  re-download per run.
- `-resume` religiously, on top of `storeDir`.
- Exclusive node allocation (`clusterOptions`) only for `bench_xlarge`.

### Compute estimate

| | Runs | CPU-hours |
|---|---|---|
| Everything except rung 3 | ~70 | ~350–450 |
| Index-only memory profile (§6) | 7 | ~5 |
| Rung 3, `quick_synteny` only | 1 | ~1 (expected OOM) |
| Rung 3, `hobrac` — optional | 1 | ~500–1,000 |

**~400 CPU-hours without rung-3 hobrac; ~1,400 with it.** Roughly 1 day
wall-clock on 144 cores for the former. RAM, not CPU, remains the binding
constraint.

---

## 8. Orchestration — a separate Nextflow pipeline

Separate repo (`synteny-bench`), not a subdirectory of `quick_synteny`. The
harness must benchmark a **pinned released version**, not the working tree,
and an independent harness is easier to defend.

```
synteny-bench/
  main.nf
  nextflow.config
  conf/{slurm.config,resources.config,cache.config}
  assets/
    datasets.csv        # accessions + md5 + truth-set pointers
    tools.csv           # tool × container digest × parameter string
  workflows/
    baseline.nf  stress_divergence.nf  stress_contiguity.nf
    stress_ploidy.nf  stress_size.nf  memory_profile.nf  simulation.nf
  modules/local/
    prepare/    fetch_assembly.nf  fetch_annotation.nf  standardize_names.nf
    cache/      busco_cached.nf  miniprot_cached.nf  blastp_cached.nf
    perturb/    fragment_assembly.nf  degrade_annotation.nf
                subset_chromosomes.nf
                simulate_rearrangements.nf  simulate_polyploid.nf
    tools/      run_quick_synteny.nf  run_hobrac.nf
                run_mcscanx.nf  run_iadhore.nf  run_jcvi.nf
                miniprot_index_only.nf        # §6 memory profile
    normalize/  normalize_<tool>.nf           # -> the §4 format
    score/      score_anchors.nf  score_chrompairs.nf
                score_homeologs.nf  score_breakpoints.nf
    collate/    collate_metrics.nf  render_report.nf
  bin/          scoring + simulation Python
```

Emits one tidy TSV — `(dataset_id, tool, tier, replicate, axis, level, subset,
metric, value, status)` — consumed by a single Quarto report. Every figure in
the paper derives from that one file.

### Resource labels — memory-tiered, not CPU-tiered

RAM is the binding constraint on this cluster, not cores.

| Label | CPUs | Memory | Use |
|---|---|---|---|
| `bench_small` | 4 | 16 GB | rung 1, scoring, perturbation |
| `bench_medium` | 12 | 48 GB | rung 2, ceiling tools |
| `bench_large` | 18 | 120 GB | memory profile, mid-size probes |
| `bench_xlarge` | 36 | 128 GB | rung 3, exclusive node |

### Failures are data

An OOM at rung 3 is a **result**, not a crash:

```groovy
errorStrategy = { task.exitStatus in [137, 140, 143] ? 'ignore' : 'retry' }
maxRetries    = 2
```

with each tool module emitting a sentinel so the normalizer records
`status=OOM` and the memory ceiling at which it occurred. The grid completes;
the failure appears as a cell in the results table.

---

## 9. The threat to the speed claim

**7.91 of hobrac's 8.46 CPU-hours is BUSCO.** Swapping in **compleasm**
(~10× faster, largely drop-in) would erase most of `quick_synteny`'s runtime
advantage.

Therefore: **benchmark hobrac both as-shipped and with compleasm**, and make
the headline comparison against the *faster* variant. If the advantage
survives, the claim is bulletproof. If it does not, better found here than in
review.

---

## 10. Reproducibility and fairness

Author-run benchmarks are suspect by default — more so with only two scored
entrants. Mitigations:

- **Default parameters for both entrants, including ours.** No tuning one side.
- **Pre-register** the protocol (OSF) before running, including §11.
- **Pin containers by digest** (`image@sha256:...`), never by tag.
- **Pin `quick_synteny` by release tag**, benchmarked as an external tool.
- Inputs declared as NCBI accessions + md5 in `assets/datasets.csv`.
- Publish `nextflow log`, all `.command.sh`, trace/report/timeline.
- Simulated truth sets + the tidy TSV deposited with a DOI.
- Fully automated, no manual steps — ideally the final grid is executed by
  someone other than the tool author.
- **Label the ceiling clearly in every figure.** It must never read as a claim
  to have beaten MCScanX or i-ADHoRe.
- Contact the hobrac authors before submission with the protocol and results.
  With n=2, getting their parameterisation confirmed is cheap insurance.

---

## 11. Pre-registered predictions

Recorded before running. Being wrong here is informative, not embarrassing.

1. Both annotation-free tools trail the MCScanX ceiling by a measurable but
   modest margin at genus/family divergence; the gap widens at order/class.
2. `quick_synteny` beats `hobrac` on anchor density and on F1 in
   poorly-BUSCO-represented phyla.
3. `quick_synteny` beats `hobrac` on runtime as-shipped; **the margin narrows
   substantially against compleasm-backed hobrac.**
4. `quick_synteny` recovers simulated homeologs and approaches the i-ADHoRe
   ceiling; `hobrac` scores ~0 by construction.
5. The MCScanX ceiling degrades sharply below ~75 % annotation completeness;
   the scored entrants are unaffected by construction.
6. `quick_synteny` OOMs at rung 3 on 128 GB; `hobrac` completes slowly.
   **Predicted loss.**
7. `quick_synteny` shows non-trivial sensitivity to proteome choice — the
   magnitude is the open question.
8. ~~The miniprot memory curve (§6) is approximately linear in genome size,
   with the 128 GB ceiling falling between 15 and 20 Gb.~~ **REFUTED** by a
   clean exclusive-node run (2026-09-20, node-c01, 384 GB, 7/7 genomes OK,
   rice through axolotl): the curve is linear (slope **9.91 GB RAM per Gb of
   genome**, forced through the origin, r fitted on 7 points 0.37–28.2 Gb) but
   the real 128 GB ceiling is **~12.9 Gb**, not 15–20. Wheat (14.6 Gb, 132 GB),
   *Pleurodeles* (19.8 Gb, 230 GB) and axolotl (28.2 Gb, 260 GB) all exceed a
   128 GB reference node despite completing cleanly here. Also notable:
   index_rss_gb ≈ align_rss_gb for every genome -- peak RAM is dominated by
   the index, not the alignment step (see APPLICATION_NOTE_PLAN.md §2).
   The FIRST attempt at this run (same genomes, shared non-exclusive nodes)
   produced nonsense across the board -- e.g. rice OOM'd at 2.82 GB against a
   16 GB request -- from node-level contention with co-scheduled jobs, not a
   real per-job limit. That data was discarded, not corrected.

---

## 12. Open decisions

- **i-ADHoRe**: painful to install (academic registration), but the only
  entrant that can contest the polyploid claim. *(Recommendation: include —
  prediction 4 is the paper's best result and needs the ceiling to mean
  anything.)*
- **Rung 3 hobrac run**: ~500–1,000 CPU-hours for one table cell, or cite
  their published scaling instead? Decide after rungs 1–2 (§7.4).
- **Chromosome subsetting (§7.2)**: acceptable for the stress axes, or does it
  weaken the result too much to be worth the saving?
- **Simulation tooling**: bespoke rearrangement simulator, or adapt an
  existing one?
- **Pre-registration venue**: OSF, or protocols.io?
