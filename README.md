# quick_synteny

Taxonomy-guided synteny plotting: given a target genome assembly and an NCBI
taxid, automatically discovers a suitable reference genome and proteome
from NCBI (climbing the taxonomy ladder from species outward until it finds
adequate assemblies), aligns the proteome against both genomes with
`miniprot`, and renders a two-genome synteny plot as a single self-contained
interactive HTML page -- a Circos-style ring, a detail panel, and a whole-
genome dotplot, all explorable and individually exportable as SVG, PNG, or JPEG.

This is a Nextflow DSL2 port of the original `legacy/quick_synteny.sh`
SLURM script, adding automatic reference-genome/proteome discovery in place
of a manually-supplied reference genome.

## Quick start

```bash
# with automatic discovery
nextflow run main.nf -profile standard \
  --taxid 562 --assembly target_genome.fa --outdir results

# with a manual reference genome/proteome (skips NCBI discovery entirely)
nextflow run main.nf -profile standard \
  --assembly target_genome.fa \
  --reference reference_genome.fa --proteome proteome.faa \
  --outdir results

# with a reference you choose on NCBI: a taxon (its best chromosome-level
# assembly) or an exact assembly; the proteome is still discovered
nextflow run main.nf -profile standard \
  --taxid 4932 --assembly target_genome.fa \
  --reference_taxid 27291 --outdir results
nextflow run main.nf -profile standard \
  --taxid 4932 --assembly target_genome.fa \
  --reference_accession GCA_056824455.1 --outdir results

# on a SLURM cluster (Apptainer/Singularity)
nextflow run main.nf -profile slurm \
  --taxid 562 --assembly target_genome.fa --outdir results
```

Needs Nextflow 25.10 or newer.

## Parameters

Every parameter is described in `nextflow_schema.json`: its type,
allowed values, default, and help. Three things read it:

- **`--help`** lists the parameters, grouped. `--help <parameter>` shows one
  parameter's full help, e.g. `nextflow run main.nf --help max_rank`.
- **Validation:** every run checks its parameters against the schema
  before anything starts. That covers types and ranges, allowed values,
  input files that must exist, and unknown parameters: a misspelled
  `--max-rank` is an error, not silently ignored. `validateParams()` in
  `main.nf` adds the rules the schema can't express: give at most one of
  `--reference`, `--reference_taxid` and `--reference_accession`, and
  `--taxid` is needed unless `--proteome` is given along with one of
  them.
- **The desktop app** takes its tips, defaults and limits from the schema.

Validation and help come from the
[nf-schema](https://github.com/nextflow-io/nf-schema) plugin, pinned in
`nextflow.config` (2.7.3, the last release for Nextflow 25.10). Nextflow
downloads it on the first run. On a machine without internet access, such
as some HPC nodes, install it once beforehand:

```bash
nextflow plugin install nf-schema@2.7.3
```

To add or change a parameter, edit both `nextflow.config`'s `params` block
and the schema. The desktop app's unit tests (`gui/test/schema.test.js`)
fail when their parameters or defaults differ. `tests/params.sh` checks
the validation and `--help`: no data needed, just Nextflow, and the
Pipeline CI workflow runs it.

## How synteny is found (no separate ortholog aligner)

The same proteome is already aligned against both genomes with `miniprot`,
inferring each protein's exon/intron structure directly on assemblies that
may have no annotation of their own. That alignment is reused directly as
the synteny signal: `miniprot` reports each protein's secondary hits
(`-N 5 --outs=0.7`, i.e. up to 5 alternative loci scoring within 70% of its
best hit) alongside its primary one, so a protein's positions in genome A
and genome B **are** the raw anchors.

`bin/extract_hits.py` turns each genome's alignments into a hit table and
groups overlapping hits into loci (hits sharing coding exons, so isoforms
count once while a locus nested inside another locus's intron stays
separate). A locus corresponds to a gene wherever the assembly's actually
annotated -- "locus" just also covers wherever the proteome aligns without
one, which is the more common case here. An **anchor** is one protein's
matched locus pair, one per genome; that's this doc's basic, rigorous unit
from here on, in place of "gene".

`bin/chain.js` then chains anchors in **rank order** rather than base
pairs, MCScanX-style: consecutive anchors of a block may skip at most
`--max_gap` anchors (default 25) on either genome, and must be strictly
increasing on both genomes (`+` blocks) or increasing on one and decreasing
on the other (`-`, inverted blocks). Counting in anchors makes the same
setting work for an anchor-dense 120 Mb genome and an anchor-sparse 30 Gb
one. A block's score is its number of anchors, i.e. distinct loci on both
genomes.

Guardrails: min identity ignores weak hits, using miniprot's Positive=
"similar-or-identical residue" score rather than its stricter Identity=.
It starts at 0.7; the page's alignment summary shows each genome's average
best-hit identity, a guide for where to move it. Min block size drops
short chains; it starts auto-tuned off the weaker genome's average identity
(15 for a close-relative pair, 0.8 or more, and 5 for a divergent one).

**Chaining is tuned after the run, on the page.** The interactive page
embeds both genomes' hit tables and re-chains in the browser (a Web Worker
running the same `bin/chain.js`), so min identity, max gap, hit rank and
min block size can all be changed live, and the blocks currently on screen
downloaded as a TSV. A run needs no chaining decisions, and the desktop
app doesn't ask for any. Every export records its settings: the TSV's
first line (`# min_identity=0.7 max_gap=25 min_block=5 ...`, named after
`bin/chain_blocks.mjs`'s flags, so it reads as the command that reproduces
the file), an SVG's `<desc>`, and every file name
(`..._blocks_id70_gap25_min5.tsv`). The pipeline still publishes
`synteny/*.links.tsv` at the starting values, with the same first line,
for use without a browser. `--min_identity`, `--max_gap` and `--min_block`
pin those starting values. They're hidden from `--help` (`--help
--showHidden` lists them, `--help min_identity` shows one), since the page
is where they're meant to be set. The exact chaining rules are specified in the header of
`bin/chain.js`.

This intentionally replaces a separate ortholog-finding aligner (an earlier
version of this pipeline used jcvi + LAST) with something cruder but much
faster and simpler: on real *Arabidopsis* genome-scale data this cut total
runtime from ~35-40 minutes to ~2-3 minutes, and removed jcvi/LAST's Apple
Silicon emulation conflict entirely (see Profiles below -- there's no longer
a Rosetta/QEMU tradeoff to make). Expect more false positives/negatives in
synteny calls than a statistically-rigorous tool like jcvi or MCScanX would
give you; this trades some of that rigor for speed and a much simpler
pipeline, which is the deliberate design goal here.

## Plot styling

Every run produces one output, `*.synteny.interactive.html` -- a single
self-contained page (Bokeh, embedded in the page -- no server or internet
connection required, open it in any browser)
with three panels side by side:

- **Ring** -- both genomes wrapped around a Circos-style circle, the
  reference genome on the top half and target on the bottom, with syntenic
  blocks drawn as ribbons crossing the middle.
- **Detail** -- click any chromosome (either genome, on the ring or on the
  dotplot's axes) to show just that chromosome's links against the other
  genome in this panel; whichever side you didn't click gets packed side by
  side. Click a single square in the dotplot's grid instead to show one
  specific (target, reference) chromosome pair, including pairs with no
  alignments at all.
- **Dotplot** -- the whole genome as a target-by-reference grid, each
  syntenic block drawn as a diagonal (or anti-diagonal, for an inversion)
  line segment.

**One selection, shown in every panel.** Click a chromosome (a ring wedge
or a dotplot axis band), a chromosome pair (a dotplot grid square) or a
single block (a ring ribbon, a dotplot segment, or a ribbon in the detail
panel), and all three panels show it: what belongs to it stays bright --
the chromosome's ribbons and partner chromosomes, its dotplot row or
column, the pair's grid square, the block itself -- while everything else
dims, and the detail panel shows the chromosome or pair. A selection
survives re-chaining, filtering and reordering for as long as it still
exists (a block that disappears falls back to its chromosome pair).
Clicking a self-link selects it on the ring only. Clear it with "Clear
selection", a click on empty ring space, or a double-click on the ring; a
click on empty detail-panel space steps a selected block back out to its
pair.

The page-wide controls sit above the panels in four groups: **Genomes**
(the two labels), **Synteny detection** (min identity, max gap, hit rank --
changing these re-chains), **Filters** (min block size, min sequence length)
and **Display** (palette, colors, chromosome order, gaps). Controls that
act on one panel alone sit under it.

A "Chromosome order" menu sets the order of the ring and both dotplot axes
alike: **Size** (the default, largest first), **File order** (each genome's
FASTA order) or **Similarity**, which places chromosomes that share blocks
next to each other -- the dotplot's synteny lines up into a clean diagonal
and the ring's ribbons untangle, most useful for a closely-related pair
with a roughly 1:1 chromosome correspondence. Both ring halves read left to
right, and the reference genome's dotplot axis reads bottom-to-top, with
the first chromosome of the chosen order at the bottom.

A "Show gaps" switch (off by default, under Display) marks assembly gaps --
runs of N's in each input FASTA, `--min_asm_gap` bp or longer (default:
100) -- on the ring and the detail panel, and as thin dotted lines on the
dotplot. Gaps are found once per genome, from the same renamed sequences
`chrom.sizes` and the links use, so their coordinates always match; each
run's gap TSVs are published to `pipeline_info/*.gaps.tsv`.

Hover any wedge, band, or ribbon for its coordinates, protein-alignment
(anchor) count, mean identity, and anchor density. Sequences are shown
under short plot names (`chr1`, `chrM`, `scaf12`, ...); an "Original
names" switch (under Display) relabels every panel, tooltip and the blocks
TSV download with each sequence's own FASTA ID instead. The full mapping is
published to `pipeline_info/*.rename_lookup.tsv`. A color-palette dropdown
and a numeric spinner recolor reference-genome chromosomes (a few curated
palettes, cycling through 1-10 discrete colors instead of one color per
chromosome); another spinner filters every panel down to blocks with at
least that many supporting anchors. Two text inputs name the two genomes
(by default their short species names, e.g. "D. melanogaster", or their
accessions when the species isn't known); the names appear in every panel
-- the titles, the ring's corner labels, the dotplot's axis titles and the
detail panel's row labels -- so any exported panel says which genome is
which. Each panel has its own save button plus a format dropdown to
export exactly that panel as SVG (a vector original -- open it in Inkscape,
Illustrator, or similar to edit it or convert it to PDF), PNG, or JPEG (both
ready to paste into a slide or document) -- named from whatever labels are
currently set. A small always-visible stats panel (under the detail panel) shows the
alignment summary -- proteome size, each genome's aligned-protein
count/mean identity, and the proteome's own origin (the species it was
auto-discovered from, or the input filename if supplied manually).

Only the **reference genome** is colored (a distinct hue per chromosome,
cycling through a qualitative palette) -- the **target** is a flat grey so
the reference genome stands out, and every link is tinted to match the
reference-genome chromosome it connects to, with opacity scaled by the
block's anchor count so strong blocks read as bolder than weak ones. This
makes it possible to trace by eye which reference-genome chromosome a
given ribbon belongs to, and to spot at a glance which blocks are
well-supported vs. marginal.

## Polyploid genomes

Every run also finds which of a genome's own chromosomes are homeologous to
each other (e.g. for an allopolyploid target or reference genome), for both
genomes, with no option to set. The page's **Show self-links** switch draws
them on the ring; it is off by default, and the page only computes them
while it is on. Self-links are drawn with a dashed outline over a fainter
fill, so they stand apart from the cross-genome synteny: a reference
genome's in its chromosomes' colors, a target genome's in dark grey. While
the switch is on, the status line next to the synteny controls counts the
self-link blocks passing Min block size (or says there are none).

No ploidy ratio needs to be declared -- the same chaining logic just runs on
each genome's hits against themselves (a protein hitting
reference-genome chr9 *and* chr10 in the same alignment already **is** the
homeolog signal), so it naturally picks up whatever multiplicity is
actually in the data. Both ends of a homeolog link land in that genome's own
half of the ring -- no special-casing needed, and the scan runs on
each genome independently (a protein forming a homeolog pair in the target
says nothing about pairing in the reference genome, and vice versa).

**Note on what the self-comparison actually shows**: on a genome with older
ancestral whole-genome-duplication history underneath the polyploidy event
you're asking about (as in Arabidopsis, which retains extensive triplicated
synteny from ancient WGDs on top of the *thaliana*/*arenosa* hybridization
that formed *suecica*), the self-scan surfaces real signal from that older
history too, not just the one hybridization event -- this is a property of
the underlying biology (and of per-anchor reciprocal matching in general), not
something this tool tries to correct for. Separating "this event's
homeologs" from "older paralogs" by age would need Ks-based dating, which is
out of scope here -- validated against a real *A. thaliana* target vs.
*A. suecica* reference genome run.

## Large genomes: chunked alignment

miniprot's whole-genome index needs about 10 GB of RAM per Gb of genome.
When that won't fit, `--miniprot_chunk_gb` splits each genome into chunks
of about that many Gb, aligns the proteome against each in parallel and
merges the results; given without a value it uses 1 Gb chunks (each task
requests `chunk Gb × --miniprot_gb_per_gb + 4` GB, ~15 GB at the default).
Whole sequences are never split, so a chunk is never smaller than its
longest chromosome.

Chunked results are **close to, but not identical to**, a whole-genome
run. miniprot decides per index which k-mers are too repetitive to seed
from and which candidate hits to keep, so a chunk sees a slightly
different picture than the whole genome does. On *C. elegans* the best hit
was identical for ~98% of proteins and ~97% of synteny blocks matched, but
~7-10% of proteins differed somewhere in their hit set -- mostly weak hits,
lower-ranked secondaries, tandem gene arrays and exact score ties -- and
the blocks that differed were all small (5-7 anchors). Each mode is
deterministic on its own.

So: leave chunking off whenever the whole-genome index fits (it's also
cheaper -- every chunk aligns the whole proteome, which cost ~4× the CPU
time on *C. elegans*); when it doesn't, use the largest chunk size that
fits, since fewer chunks stay closer to whole-genome; and compare runs only
at the same setting.

## Profiles

- `standard` -- Docker, local executor (laptop/CI).
- `conda` -- no container engine: conda environments built from
  `envs/*.yml`, local executor. Needs conda, mamba or micromamba (add
  `-c` config with `conda.useMicromamba = true` for the latter).
- `slurm` -- Apptainer/Singularity, SLURM executor (HPC).
- `test` -- Docker, local executor, tiny resource caps.

`conda` pins the same tool versions as the images, but each host solves
those environments itself, so dependencies and build strings can differ
from the images; results can too, though on an *S. cerevisiae* run every
output matched the `standard` profile's byte for byte. Compare runs made
with the same profile.

In the container profiles, every process runs a Seqera Containers image
for the host's architecture (linux/arm64 on Apple Silicon, linux/amd64
otherwise) -- no `--platform` is forced and nothing runs emulated. That
matters beyond speed: node's Maglev compiler miscompiles the chainer under
x86-64 emulation, silently skewing its auto-tuned parameters.
`bin/chain_blocks.mjs` also disables Maglev and fails (exit 3, retried) on
an impossible result, in case it ever does run emulated.

## Containers

Tool versions live in one place: `envs/*.yml`, one conda environment per
image (ncbi-datasets + unzip, python, miniprot, nodejs, bokeh). Each process
names its environment with an `env_<name>` label. `conf/containers.config`
maps each label to its four images -- Docker and Singularity (an HTTPS SIF
download), each for linux/amd64 and linux/arm64 -- and is generated, not
edited by hand:

```bash
tools/update_containers.py            # every environment
tools/update_containers.py miniprot   # just envs/miniprot.yml
```

It asks the public Wave API behind seqera.io/containers (no token needed)
for each image, waits for any that have to be built first, and rewrites
the config. The same package list always gives the same frozen image, so a
rerun only builds what changed. To update a tool, bump its version in
`envs/<name>.yml`, run the script, and commit both files. Runs themselves
never contact Wave -- the image URLs are fixed in the repo.

## Repository layout

- `main.nf`, `nextflow.config`, `conf/` -- pipeline entrypoint and profiles.
- `nextflow_schema.json` -- every parameter: validation, `--help`, and the
  desktop app's form.
- `tests/params.sh` -- parameter validation and `--help` checks.
- `modules/local/` -- one process (or small process group) per tool.
- `bin/` -- Python/shell helper scripts used inside processes.
- `envs/`, `tools/update_containers.py` -- container environments and the
  script that generates `conf/containers.config` from them.
- `legacy/` -- the original bash/SLURM script this pipeline replaces.
- `gui/` -- the desktop app (see `gui/README.md`).
