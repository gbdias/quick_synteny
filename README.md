# quick_synteny

Taxonomy-guided synteny plotting: given a target genome assembly and an NCBI
taxid, automatically discovers a suitable reference genome and proteome
from NCBI (climbing the taxonomy ladder from species outward until it finds
adequate assemblies), aligns the proteome against both genomes with
`miniprot`, and renders a two-genome synteny plot as a single self-contained
interactive HTML page -- a Circos-style ring, a zoom panel, and a whole-
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

# on a SLURM cluster (Apptainer/Singularity)
nextflow run main.nf -profile slurm \
  --taxid 562 --assembly target_genome.fa --outdir results
```

See `nextflow run main.nf --help` for the full parameter list.

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

Guardrails: `--min_identity` ignores weak hits. By default it's auto-tuned
to the weaker genome's mean best-hit identity, `max(0.3, min(0.9,
weaker_mean))`, using miniprot's Positive= "similar-or-identical residue"
score rather than its stricter Identity=, so a divergent pair isn't forced
through a threshold tuned for close relatives. `--min_block` drops short
chains; it's auto-tuned off the same number (15 for a close-relative pair,
5 for a divergent one). Pass either to pin it yourself. These only set the
starting point: the interactive page embeds both genomes' hit tables and
re-chains in the browser (a Web Worker running the same `bin/chain.js`), so
min identity, max gap, hit rank and min block size can all be changed live,
and the blocks currently on screen downloaded as a TSV. The exact chaining
rules are specified in the header of `bin/chain.js`.

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
self-contained page (Bokeh, no server required -- open it in any browser)
with three panels side by side:

- **Ring** -- both genomes wrapped around a Circos-style circle, the
  reference genome on the top half and target on the bottom, with syntenic
  blocks drawn as ribbons crossing the middle.
- **Zoom** -- click any chromosome (either genome, on the ring or on the
  dotplot's axes) to zoom this panel into just that chromosome's links
  against the other genome; whichever side you didn't click gets packed
  side by side. Click a single square in the dotplot's grid instead to zoom
  straight into one specific (target, reference) chromosome pair, including
  pairs with no alignments at all.
- **Dotplot** -- the whole genome as a target-by-reference grid, each
  syntenic block drawn as a diagonal (or anti-diagonal, for an inversion)
  line segment. An "order chromosomes by similarity" toggle reorders both
  axes so shared synteny lines up into a clean diagonal -- most useful for a
  closely-related pair with a roughly 1:1 chromosome correspondence.

An "Order by size" switch (on by default) sorts every chromosome -- on the
ring and both dotplot axes alike -- largest to smallest; switch it off to
see each genome's chromosomes in their original FASTA order instead. The
reference genome's dotplot axis always reads bottom-to-top with the
largest (or, with the switch off, the first-in-file) chromosome at the
bottom, under either setting.

A "Show gaps" switch (off by default, below the ring) marks assembly gaps --
runs of N's in each input FASTA, `--min_asm_gap` bp or longer (default:
100) -- on the ring and the zoom panel, and as thin dotted lines on the
dotplot. Gaps are found once per genome, from the same renamed sequences
`chrom.sizes` and the links use, so their coordinates always match; each
run's gap TSVs are published to `pipeline_info/*.gaps.tsv`.

Hover any wedge, band, or ribbon for its coordinates, protein-alignment
(anchor) count, mean identity, and anchor density. A color-palette dropdown
and a numeric spinner recolor reference-genome chromosomes (a few curated
palettes, cycling through 1-10 discrete colors instead of one color per
chromosome); another spinner filters every panel down to blocks with at
least that many supporting anchors. Two text inputs relabel "target"/
"reference" to the actual species/genome names everywhere a title or bar
shows them, and each panel has its own save button plus a format dropdown to
export exactly that panel as SVG (a vector original -- open it in Inkscape,
Illustrator, or similar to edit it or convert it to PDF), PNG, or JPEG (both
ready to paste into a slide or document) -- named from whatever labels are
currently set. A small always-visible stats panel (under the zoom panel) shows the
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
while it is on.

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

## Profiles

- `standard` -- Docker, local executor (laptop/CI).
- `slurm` -- Apptainer/Singularity, SLURM executor (HPC).
- `test` -- Docker, local executor, tiny resource caps.

All containers now run fine under Docker Desktop's default Rosetta
emulation on Apple Silicon -- an earlier jcvi/LAST-based version of this
pipeline needed a QEMU/Rosetta tradeoff that no longer applies now that
jcvi has been removed from the pipeline entirely.

## Repository layout

- `main.nf`, `nextflow.config`, `conf/` -- pipeline entrypoint and profiles.
- `modules/local/` -- one process (or small process group) per tool.
- `bin/` -- Python/shell helper scripts used inside processes.
- `legacy/` -- the original bash/SLURM script this pipeline replaces.
