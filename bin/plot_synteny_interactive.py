#!/usr/bin/env python3
"""Render a two-genome synteny plot as a self-contained interactive HTML page
(Bokeh, no server -- all interactivity is CustomJS). This is the pipeline's
only plot output -- there is deliberately no static-image equivalent, so
that "explore, then export exactly the panel you want" (see each panel's
save button + format selector, SYN.exportFigure -- SVG for a vector original
to edit or convert to PDF, PNG for something ready to paste into a slide or
document) can be the one workflow, rather than needing separate static
renders kept in sync with it.

Three-column layout: a full Circos-style ring on the left (every wedge/
ribbon individually hoverable), a linear zoom/detail panel in the middle, a
whole-genome dotplot on the right. Clicking a reference-chromosome wedge, a
target-chromosome wedge, or a chromosome band on either dotplot axis zooms
the middle panel into that one chromosome's links -- whichever side you
didn't click gets packed (multiple chromosomes side by side if more than one
is involved), the clicked chromosome stays alone on its own row. Clicking a
single square in the dotplot grid instead zooms straight into that one
(target, reference) chromosome pair -- always exactly two bars, no packing,
and this is also the only way to select a pair with zero links (see
SYN.buildPairDetailData), since the axis rulers can only ever select one
chromosome at a time. Reference-genome chromosomes are always drawn on top
in every panel, regardless of which side triggered the zoom. A dropdown
picks which of a few curated color palettes (see PALETTES) reference-genome
chromosomes cycle through, and a numeric spinner recolors every panel by
cycling through 1-MAX_COLORS discrete colors from whichever one is active
instead of one color per chromosome; a min-block-size spinner filters every
panel down to blocks with at least that many supporting anchors -- since
extraction is independent of the minimum block size (docs/specs/hit_table.md
section 4.6), the chainer is always run at minBlock=3 and this spinner is a
pure client-side filter on the result (see SYN.data.ribbons, built by
SYN.buildRingLayout from whatever SYN.applyChainResult most recently put in
SYN.data.linksByQuery/linksByReference, and SYN.buildOverviewRibbons). Min
identity, max gap, and hit rank each trigger a fresh client-side re-chain
(SYN.requestChain) -- see SYN.startChainer/SYN.applyChainResult below. Spinners
rather than sliders since the filters are plain numeric comparisons that work
for any value, not just a handful of steps -- a free-form numeric input
doesn't imply a false ceiling the way a slider's end-of-track does. Two text
inputs let a viewer relabel "target"/"reference" to the actual genome/species
names before exporting a panel (see SYN.applyLabels) -- every exportable
title/bar-label reads from these rather than the pipeline's literal role
tags. A stats panel (see SYN.applyStats) shows the alignment summary --stats
optionally provides, since that no longer has anywhere to live inside an
exported image the way it did on the static plots this replaced.

Every wedge/ribbon/segment/gap polygon on the page -- ring, zoom panel, and
dotplot alike -- is built once, client-side, as plain JS (SHARED_JS): SYN.init
(see build_page()'s doc.js_on_event(DocumentReady, ...)) runs it for the
initial render, through the exact same functions every later control change
reuses, so the browser can rebuild an arbitrary (pivot side, pivot chromosome,
color count, min block size, chromosome order) combination instead of only
combinations precomputed ahead of time. Python never computes any geometry,
nor any synteny block, at all: it only ships each genome's raw hit table
(docs/specs/hit_table.md section 6, SYN.hitsPayload) plus chromosome sizes
and palette indices, as embedded JSON, and inlines bin/chain.js so the browser
can chain (and re-chain, on every parameter change) client-side, in a Web
Worker when one is available -- see SYN.startChainer/SYN.applyChainResult.

This is not run through bin/'s usual container -- it needs Bokeh, which has
no bioconda recipe, so this runs on a Seqera Containers (Wave) image built
directly from conda-forge's bokeh package instead. See
modules/local/pygenomeviz_plot.nf for the full reasoning.
"""
import argparse
import base64
import gzip
import html
import json
import math
import os
import sys

import numpy as np

from bokeh.document import Document
from bokeh.embed import file_html
from bokeh.events import DocumentReady, DoubleTap, Tap
from bokeh.layouts import column, row
from bokeh.models import (ColumnDataSource, HoverTool, TapTool, CustomJS, CustomJSTickFormatter,
                           Button, Div, Range1d, Select, Spinner, Switch, TextInput, Toggle)
from bokeh.plotting import figure
from bokeh.resources import CDN

PALETTE = [
    '#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B2',
    '#937860', '#DA8BC3', '#8C8C8C', '#CCB974', '#64B5CD',
]
# Okabe-Ito (the first 8) plus two extra hues, for a colorblind-friendlier
# alternative to PALETTE's seaborn "deep" set above
PALETTE_COLORBLIND = [
    '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
    '#D55E00', '#CC79A7', '#000000', '#6A3D9A', '#B15928',
]
# seaborn's "pastel" set -- same hues as PALETTE, softened
PALETTE_PASTEL = [
    '#A1C9F4', '#FFB482', '#8DE5A1', '#FF9F9B', '#D0BBFF',
    '#DEBB9B', '#FAB0E4', '#CFCFCF', '#FFFEA3', '#B9F2F0',
]
# seaborn's "bright" set -- same hues as PALETTE, more saturated/vivid
PALETTE_BRIGHT = [
    '#023EFF', '#FF7C00', '#1AC938', '#E8000B', '#8B2BE2',
    '#9F4800', '#F14CC1', '#A3A3A3', '#FFC400', '#00D7FF',
]
# seaborn's "dark" set -- same hues again, deepened/desaturated; pairs well
# with light text or a light background where PALETTE/PALETTE_BRIGHT read
# as too loud
PALETTE_DARK = [
    '#001C7F', '#B1400D', '#12711C', '#8C0800', '#591E71',
    '#592F0D', '#A23582', '#3C3C3C', '#B8850A', '#006374',
]
# matplotlib's "tab10" -- a different, widely-recognized categorical set
# (distinct hues from seaborn's "deep"/PALETTE above, not just a
# lighter/darker variant of the same ten) for a viewer who wants a familiar
# look from outside the seaborn family
PALETTE_TABLEAU = [
    '#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD',
    '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF',
]
# every option needs exactly MAX_COLORS entries -- SYN.paletteColor (see
# SHARED_JS) indexes with idx % k for any k up to MAX_COLORS, so a shorter
# list would read past its own end and recolor as undefined
PALETTES = {
    'Default': PALETTE,
    'Colorblind-safe': PALETTE_COLORBLIND,
    'Pastel': PALETTE_PASTEL,
    'Bright': PALETTE_BRIGHT,
    'Dark': PALETTE_DARK,
    'Tableau': PALETTE_TABLEAU,
}
DEFAULT_PALETTE_NAME = 'Default'
TARGET_GREY = '#999999'

OUTER_R = 1.0
RING_WIDTH = 0.045
INNER_R = OUTER_R - RING_WIDTH
LINK_R = INNER_R
GROUP_GAP = math.radians(4)
CHROM_GAP = math.radians(0.5)

# an assembly gap is typically a few hundred bp against a multi-Mb
# chromosome -- drawn at its true angular width, most gaps would render
# sub-pixel and be both invisible and practically unhoverable. Widening only
# the ones that would otherwise fall under this floor (real length is still
# exact in the tooltip/label) is the same "minimum feature width" convention
# other genome browsers use for tiny features -- not applied to a gap that's
# already wide enough to show its true size. Small on purpose: a gap wedge
# spans the ideogram's FULL radial height (INNER_R to OUTER_R, same as
# SYN.wedgePolygonJS's own r0/r1 defaults -- no separate radial inset), so it should
# read as a thin needle mark across that whole band, not a short wide block --
# only the angular width (this constant) needs flooring, not the height.
MIN_GAP_ANGLE = math.radians(0.2)

# the zoom panel figure's own fixed pixel width -- embedded (not duplicated
# as a separate JS literal) so SYN.widenGapBp's px-to-bp conversion can
# never drift from the actual figure build_page() creates below
DETAIL_FIG_WIDTH = 446

# the default right-hand toolbar Bokeh reserves on every figure here (pan/
# wheel-zoom/reset/hover icons) -- that strip eats into detail_fig's own
# DETAIL_FIG_WIDTH rather than sitting outside it, so the panel's actual
# *drawn* frame (what a viewer reads as "the zoom panel") is this much
# narrower on the right than the figure's nominal width.
PLOT_TOOLBAR_WIDTH = 30

# the stats box below detail_fig -- matched to the zoom panel's drawn frame,
# not its full nominal width, so the box's own right edge lines up with the
# frame's right edge instead of running on past it under the toolbar icons
STATS_PANEL_WIDTH = DETAIL_FIG_WIDTH - PLOT_TOOLBAR_WIDTH

# every titleless toolbar control below a figure (save/clear buttons, the
# SVG/PNG/JPEG format dropdowns, the "order by similarity" toggle) -- left
# unset, each one auto-sizes to its own label/font metrics instead (a plain
# Button came out a couple px shorter than one with an emoji glyph in its
# label, and Select's own browser-native control renders shorter still than
# either), so the buttons along one row didn't line up. Explicit and shared
# so every one of them matches, no matter what's in its label.
TOOLBAR_CONTROL_HEIGHT = 32

# the header row's titled controls (target/reference label, color palette,
# colors, min identity, max gap, hit rank, min block size) -- narrower than
# their old 220px, which was much wider than any of their labels or values
# actually need
TOP_CONTROL_WIDTH = 150

# dotplot gap lines only -- a lighter grey than the ring/zoom panel's own
# solid-black gap markers (a `multi_line` renderer style, set once here in
# Python; no JS mirror needed since it's never re-set per-point client-side)
GAP_LINE_COLOR = '#AAAAAA'

MAX_COLORS = 10
DEFAULT_COLORS = MAX_COLORS

# fraction of each axis's total span given to the dotplot's chromosome-ruler
# strips (see SHARED_JS's SYN.buildDotplotLayout) -- the clickable regions
# standing in for x/y axis ticks, since Bokeh has no native "clickable tick
# label"
DP_RULER_FRAC = 0.035

def read_chrom_sizes(path):
    chroms = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            name, size = line.split('\t')[:2]
            chroms.append((name, int(size)))
    return chroms


def read_gaps(path):
    """bin/rename_sequences.py's --out_gaps output: chrom, start, end -- no
    header, 0-based half-open. Returns [] for path=None (no --target_gaps/
    --comparison_gaps given)."""
    if not path:
        return []
    gaps = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            chrom, start, end = line.split('\t')[:3]
            gaps.append({'chrom': chrom, 'start': int(start), 'end': int(end)})
    return gaps


HITS_TSV_HEADER = ['chrom', 'start', 'end', 'strand', 'positive', 'identity', 'score', 'rank', 'protein', 'locus']


def read_hits_tsv(raw_bytes):
    """Parses one bin/extract_hits.py output (docs/specs/hit_table.md section
    1) already read into memory as raw gzip bytes -- into per-column lists,
    in TSV row order (already sorted by (chrom, start, end, protein, rank),
    see that script). `identity`/`score` aren't part of the in-memory
    HitTable (spec section 2) and are dropped here."""
    text = gzip.decompress(raw_bytes).decode('utf-8')
    lines = text.split('\n')
    header = lines[0].split('\t')
    if header != HITS_TSV_HEADER:
        sys.exit(f"ERROR: unexpected hits TSV header {header}")
    chrom, start, end, strand, positive, rank, protein, locus = [], [], [], [], [], [], [], []
    for line in lines[1:]:
        if not line:
            continue
        row = line.split('\t')
        chrom.append(row[0]); start.append(int(row[1])); end.append(int(row[2]))
        strand.append(row[3]); positive.append(float(row[4])); rank.append(int(row[7]))
        protein.append(row[8]); locus.append(int(row[9]))
    return {'chrom': chrom, 'start': start, 'end': end, 'strand': strand,
            'positive': positive, 'rank': rank, 'protein': protein, 'locus': locus}


def encode_hits_column(values, dtype):
    """One GenomeColumns column (docs/specs/hit_table.md section 6): little-
    endian typed array, gzip'd, base64'd. numpy is used only for the
    explicit little-endian dtypes and fast bulk casting -- it ships in the
    bokeh container (see modules/local/pygenomeviz_plot.nf)."""
    np_dtype = {'u8': '<u1', 'u16': '<u2', 'u32': '<u4', 'f64': '<f8'}[dtype]
    data = np.asarray(values, dtype=np_dtype).tobytes()
    return {'dtype': dtype, 'data': base64.b64encode(gzip.compress(data, mtime=0)).decode('ascii')}


def widened_dtype(values, narrow, wide, limit):
    """spec section 6: startDelta/length are u32, f64 "if any value >= 2**32"
    (a chromosome bigger than 4.29 Gb)."""
    return wide if values and max(values) >= limit else narrow


def build_genome_columns(hits, chrom_names, prot_index):
    """One GenomeColumns dict (docs/specs/hit_table.md section 6) for a
    single genome's parsed hit table (read_hits_tsv). chrom_names is that
    genome's FULL chrom.sizes-order name list (Dataset's own *_chroms_natural,
    not filtered by --min_seq_size or anything else) -- read_hits_tsv's rows
    were already restricted to sequences in that same chrom.sizes by
    bin/extract_hits.py, so every chrom here is guaranteed present."""
    chrom_index = {name: i for i, name in enumerate(chrom_names)}
    n = len(hits['chrom'])
    chrom_idx = [chrom_index[c] for c in hits['chrom']]

    # startDelta (spec section 6): start minus the previous row's start on
    # the same chrom, first row per chrom relative to 0 -- rows are already
    # grouped by chrom (TSV row order, see read_hits_tsv's docstring)
    start_delta = [0] * n
    prev_chrom, prev_start = None, 0
    for i in range(n):
        if hits['chrom'][i] != prev_chrom:
            prev_chrom, prev_start = hits['chrom'][i], 0
        start_delta[i] = hits['start'][i] - prev_start
        prev_start = hits['start'][i]
    length = [hits['end'][i] - hits['start'][i] for i in range(n)]
    strand = [1 if s == '+' else 0 for s in hits['strand']]
    positive = [round(p * 10000) for p in hits['positive']]
    rank = [min(r, 255) for r in hits['rank']]
    locus = hits['locus']
    n_loci = (max(locus) + 1) if locus else 0

    return {
        'n': n, 'nLoci': n_loci,
        'cols': {
            'chrom': encode_hits_column(chrom_idx, 'u16'),
            'startDelta': encode_hits_column(start_delta, widened_dtype(start_delta, 'u32', 'f64', 2 ** 32)),
            'length': encode_hits_column(length, widened_dtype(length, 'u32', 'f64', 2 ** 32)),
            'strand': encode_hits_column(strand, 'u8'),
            'positive': encode_hits_column(positive, 'u16'),
            'rank': encode_hits_column(rank, 'u8'),
            'prot': encode_hits_column([prot_index[p] for p in hits['protein']], 'u32'),
            'locus': encode_hits_column(locus, 'u32'),
        },
    }


def build_hits_payload(target_hits_path, comparison_hits_path, target_chrom_names, comparison_chrom_names):
    """SYN.hitsPayload (docs/specs/hit_table.md section 6) -- the one thing
    Python ships instead of precomputed links/blocks. `reference:
    'same_as_target'` (skipping a second, redundant copy of the same bytes)
    only when both the raw hit tables AND the two genomes' chrom.sizes
    orderings are byte-for-byte identical -- SYNCHAIN.decodePayload's
    same_as_target path reuses the target table's chrom-index column
    verbatim under the reference's own chromNames/chromSizes, which is only
    correct if chromosome N means the same thing on both sides."""
    with open(target_hits_path, 'rb') as f:
        target_bytes = f.read()
    same_genome = target_chrom_names == comparison_chrom_names
    if same_genome:
        with open(comparison_hits_path, 'rb') as f:
            same_genome = f.read() == target_bytes

    target_hits = read_hits_tsv(target_bytes)
    comparison_hits = None if same_genome else read_hits_tsv(open(comparison_hits_path, 'rb').read())

    # union of both genomes' accessions, sorted lexicographically (spec
    # section 2) -- ASCII protein accessions sort identically under Python's
    # code-point order and JS's default (UTF-16 code unit) Array.sort()
    all_proteins = set(target_hits['protein']) | (set(comparison_hits['protein']) if comparison_hits else set())
    proteins = sorted(all_proteins)
    prot_index = {p: i for i, p in enumerate(proteins)}

    target_cols = build_genome_columns(target_hits, target_chrom_names, prot_index)
    reference_cols = ('same_as_target' if same_genome
                       else build_genome_columns(comparison_hits, comparison_chrom_names, prot_index))

    proteins_gz = gzip.compress('\n'.join(proteins).encode('utf-8'), mtime=0)
    return {
        'version': 1,
        'proteins': base64.b64encode(proteins_gz).decode('ascii'),
        'genomes': {'target': target_cols, 'reference': reference_cols},
    }


class Dataset:
    def __init__(self, query_chrom_sizes, subject_chrom_sizes, query_gaps_path=None, subject_gaps_path=None):
        # RENAME_SEQUENCES emits natural (FASTA) order now, not size order
        # (see main.nf's own comment) -- self.query_chroms/subject_chroms
        # stay the name used everywhere below (dotplot layout, color
        # assignment) and represent whichever order is DEFAULT-active, which
        # is size order (size_order_toggle defaults to on in build_page() --
        # kept in sync deliberately, same as every other Python/JS
        # default-state pair in this file). The natural-order lists are kept
        # alongside purely so build_page() can embed them for the
        # client-side reorder engines (SYN.buildRingLayout, SYN.dpNaturalOrder)
        # to switch back to when that toggle goes off, AND so SYN.startChainer
        # has the full (unfiltered) chromosome name/size lists the embedded
        # hit tables' chrom indices refer to (see build_hits_payload).
        self.query_chroms_natural = read_chrom_sizes(query_chrom_sizes)
        self.subject_chroms_natural = read_chrom_sizes(subject_chrom_sizes)
        self.query_chroms = sorted(self.query_chroms_natural, key=lambda c: -c[1])
        self.subject_chroms = sorted(self.subject_chroms_natural, key=lambda c: -c[1])
        self.query_sizes = dict(self.query_chroms)
        self.subject_sizes = dict(self.subject_chroms)

        # assembly gaps -- filtered to chroms actually present, same
        # reasoning links used to get filtered by before blocks moved
        # client-side
        self.query_gaps = [g for g in read_gaps(query_gaps_path) if g['chrom'] in self.query_sizes]
        self.subject_gaps = [g for g in read_gaps(subject_gaps_path) if g['chrom'] in self.subject_sizes]


def group_gaps_by_chrom(gaps):
    """{chrom: [{start, end}, ...]} -- the form the ring, the zoom panel, and
    the dotplot's gap lines all need (per-chromosome lookup, not a flat
    wedge-polygon list), embedded once and shared by all three client-side
    (see SHARED_JS's targetGapsByChrom/referenceGapsByChrom usage, read by
    SYN.buildGapRecords/SYN.gapsForChrom/SYN.buildDotplotGapLines)."""
    by_chrom = {}
    for g in gaps:
        by_chrom.setdefault(g['chrom'], []).append({'start': g['start'], 'end': g['end']})
    return by_chrom


def trunc1(x):
    """Truncates (not rounds) to 1 decimal place -- e.g. 99.96 -> 99.9, not
    100.0, so the stats panel never implies a percentage reached a round
    number it didn't actually reach. Mirrors SHARED_JS's own truncation in
    SYN.buildStatsHtml (Math.floor(x*10)/10) -- kept in sync deliberately."""
    return math.floor(x * 10) / 10


def format_stats_html(target_label, reference_label, stats):
    """The interactive HTML's replacement for the (removed) static plots'
    stats inset -- see compute_alignment_stats.py. Mirrors SHARED_JS's
    SYN.buildStatsHtml exactly (same markup), since that's what re-renders
    this same panel client-side whenever a viewer edits the target/reference
    labels -- kept in sync deliberately (see module docstring)."""
    if not stats:
        return ''
    total = stats['proteome_total']
    origin = stats.get('proteome_origin')

    def row(label, aligned, identity):
        pct = (aligned / total * 100) if total else 0.0
        return (f"<tr><td style='padding-right:20px'>{html.escape(label)}</td>"
                f"<td style='text-align:right;padding-right:20px'>{aligned:,} ({trunc1(pct):.1f}%)</td>"
                f"<td style='text-align:right'>{trunc1(identity * 100):.1f}%</td></tr>")

    origin_line = (f"<div style='color:#666;margin-top:2px'>Proteome: {html.escape(origin)}</div>"
                   if origin else '')

    # display:block + an explicit px width (matching stats_div's own
    # width=STATS_PANEL_WIDTH, box-sizing:border-box so the border/padding
    # don't push it wider) -- rather than the old shrink-to-fit
    # display:inline-block -- makes this box span the zoom panel's drawn
    # frame instead of only as wide as its own text/table content.
    return (
        "<div style='font-size:13px;color:#333;border:1px solid #ddd;border-radius:6px;"
        f"padding:8px 12px;display:block;box-sizing:border-box;width:{STATS_PANEL_WIDTH}px;"
        "margin:4px 0 8px 0'>"
        f"<b>Alignment summary</b><br>{total:,} input proteins"
        f"{origin_line}"
        "<table style='border-collapse:collapse;margin-top:4px'>"
        "<tr style='color:#666'><th style='text-align:left;padding-right:20px'></th>"
        "<th style='text-align:right;padding-right:20px'>aligned</th>"
        "<th style='text-align:right'>avg identity</th></tr>"
        + row(target_label, stats['query_aligned'], stats['query_mean_identity'])
        + row(reference_label, stats['subject_aligned'], stats['subject_mean_identity'])
        + "</table></div>"
    )


def build_overview_sources():
    """The ring's four ColumnDataSources -- query/subject wedges, ribbons,
    and chromosome labels -- created with their real columns but no rows.
    Python computes no wedge/ribbon geometry at all any more (see module
    docstring): SYN.init (SHARED_JS, run from build_page()'s
    doc.js_on_event(DocumentReady, ...)) fills all four via
    SYN.applyRingLayout as soon as the page's document is ready, through the
    exact same function every later reorder/filter/recolor reuses."""
    q_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[], name=[], size_label=[], group=[]))
    s_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[], name=[], size_label=[], group=[], palette_index=[]))
    r_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[], alpha=[], label=[], palette_index=[]))
    label_src = ColumnDataSource(dict(x=[], y=[], text=[]))
    return q_src, s_src, r_src, label_src


def build_dotplot_sources():
    """The dotplot's order-dependent ColumnDataSources -- query/subject
    ruler strips, block segments, boundary gridlines, and axis labels --
    created with their real columns but no rows, same reasoning as
    build_overview_sources() above. SYN.init fills all six via
    SYN.applyDotplotOrder once the document is ready, which also sets the
    figure's actual x_range/y_range (see build_page()'s placeholder
    Range1ds) -- there is no longer an invisible grid-cell hit layer here at
    all (see SYN.dpCellAt, which hit-tests by binary search instead)."""
    dp_q_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[], name=[], size_label=[], group=[]))
    dp_s_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[], name=[], size_label=[], group=[], palette_index=[]))
    dp_seg_src = ColumnDataSource(dict(xs=[], ys=[], line_color=[], alpha=[], label=[], palette_index=[]))
    dp_grid_src = ColumnDataSource(dict(xs=[], ys=[]))
    dp_q_label_src = ColumnDataSource(dict(x=[], y=[], text=[]))
    dp_s_label_src = ColumnDataSource(dict(x=[], y=[], text=[]))
    return dp_q_src, dp_s_src, dp_seg_src, dp_grid_src, dp_q_label_src, dp_s_label_src


# JS port of the linear-layout math above (an angle-free version of the ring's
# angular offsets + the bar/ribbon polygon construction), generalized so
# either side can be the "pivot" (the chromosome that was clicked, drawn
# alone on its own row) with the other side packed -- see module docstring
# for why this is the only implementation of any of it any more.
SHARED_JS = r"""
window.SYN = window.SYN || {};
// mode is 'pivot' (a ring wedge or dotplot ruler was clicked -- pivotSide/
// pivotName apply) or 'pair' (a dotplot grid cell was clicked -- pairTarget/
// pairSubject apply). The color/min-block-size spinners need to know which
// of SYN.buildDetailData/SYN.buildPairDetailData to re-run.
// dpQueryOffsets/dpSubjectOffsets are the dotplot's *current* per-chromosome
// axis offsets (natural order at load, replaced wholesale by the "order by
// similarity" toggle -- see SYN.applyDotplotOrder) -- every dotplot redraw
// that isn't a reorder itself (recoloring, the min-block-size filter) reads
// these rather than recomputing an order, so it always draws whatever order
// is currently active without needing to know which one that is.
// targetLabel/referenceLabel are the viewer-editable display names (see
// SYN.applyLabels) -- every title/bar-label that would otherwise show the
// pipeline's literal "target"/"reference" role tags reads from these
// instead, seeded to those same literal tags at load (see SYN.init) so
// nothing changes on screen until a viewer actually types something into
// the label inputs.
SYN.state = {
    pivotSide: null, pivotName: null, mode: null, pairTarget: null, pairSubject: null,
    dpQueryOffsets: null, dpSubjectOffsets: null,
    // the dotplot's current per-axis chromosome order itself (not just the
    // offsets above) -- SYN.dpCellAt binary-searches these to hit-test a tap
    // without a per-cell renderer (see SYN.applyDotplotOrder, the only
    // place either is ever set)
    dpQueryOrder: null, dpSubjectOrder: null,
    targetLabel: null, referenceLabel: null,
    // the zoom panel's own x_range/y_range are rewritten on every click (see
    // SYN.applyDetail below), unlike the ring/dotplot's fixed ranges -- so
    // its double-click-to-reset target has to track whatever was most
    // recently set programmatically, not a single fixed value baked in at
    // page load (see build_page()'s three DoubleTap handlers)
    detailRange: null,
    // a snapshot of the ring ribbons' alpha array taken the moment a ribbon
    // is first clicked, so SYN.applyRibbonHighlight can restore the exact
    // pre-highlight (score-based) alphas on deselect rather than guessing --
    // null whenever no ribbon is currently highlighted
    ribbonBaseAlpha: null,
    // read by SYN.dpNaturalOrder (both the dotplot's own toggle and the
    // ring's reorder read through that one function) -- true = size order
    // (SYN.data.queryNames/subjectNames), false = FASTA/natural order
    // (SYN.data.queryNamesNatural/subjectNamesNatural). Matches
    // size_order_toggle's default (see build_page()) -- kept in sync
    // deliberately, same as every other default-state pair in this file.
    orderBySize: true,
    // read by SYN.applyChainResult, to decide whether a fresh chain result
    // needs a full dotplot reorder (similarity order depends on link
    // weights, which just changed) or only a segment redraw in place.
    // Matches order_toggle's default (see build_page()) -- kept in sync
    // deliberately, same as every other default-state pair in this file.
    orderBySimilarity: false,
    // read by SYN.applyDetail (gates whether a freshly built detail's gap
    // ticks actually reach detail_gap_source) and by show_gaps_toggle's own
    // callback (ring + dotplot gap visibility) -- matches show_gaps_toggle's
    // default (see build_page()) -- kept in sync deliberately, same as
    // every other default-state pair in this file.
    showGaps: false,
    // bp floor read by SYN.filterBySize -- chromosomes shorter than this are
    // dropped from every panel (see min_seq_size_spinner in build_page()).
    // 0 by default: --min_seq_size already dropped anything smaller than the
    // pipeline's own threshold before this script ever saw the data (see
    // RENAME_SEQUENCES's chrom_sizes output, main.nf), so 0 here means "show
    // everything actually embedded in this page", not "no filtering ever
    // happened upstream".
    minSeqSize: 0,
    // the dotplot's current total span and ruler-strip thickness -- the
    // min-length filter above can shrink these (fewer/smaller chromosomes ->
    // smaller total) even for the same order, so SYN.buildDotplotLayout
    // recomputes them on every call and SYN.applyDotplotOrder refreshes
    // these four alongside the offsets. Set for real by SYN.init's own
    // first call to SYN.applyDotplotOrder, before the dotplot's
    // DoubleTap-to-reset handler or SYN.buildDotplotGapLines can ever run.
    dpTotalX: null, dpTotalY: null, dpRx: null, dpRy: null,
};

// how dim the OTHER ribbons go while one is highlighted -- low enough to
// read as "not what you clicked", not so low they vanish and lose all
// spatial context for where the highlighted one sits among the rest
SYN.RIBBON_DIM_ALPHA = 0.08;

// The label TextInputs are free-text (see build_page()'s target_label_input/
// reference_label_input) and, unlike chromosome names, are never validated
// against anything -- every place their value gets interpolated into HTML
// rendered via a Div's .text or a HoverTool's "{safe}" formatter (both raw
// HTML, not escaped by Bokeh) runs it through this first, so a label
// containing e.g. "<" can't break the page's markup.
SYN.escapeHtml = function(s) {
    return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;').replace(/'/g, '&#39;');
};

SYN.lighten = function(hex, amount) {
    amount = amount === undefined ? 0.35 : amount;
    const r = parseInt(hex.slice(1, 3), 16), g = parseInt(hex.slice(3, 5), 16), b = parseInt(hex.slice(5, 7), 16);
    const blend = (c) => Math.round(c + (255 - c) * amount);
    const toHex = (c) => c.toString(16).padStart(2, '0');
    return '#' + toHex(blend(r)) + toHex(blend(g)) + toHex(blend(b));
};

SYN.paletteColor = function(index, k) {
    return SYN.data.palette[index % k];
};

// black or white, whichever contrasts better against a given ideogram fill
// color -- needed because some palettes (e.g. Colorblind-safe/Okabe-Ito)
// include near-black swatches, and the zoom panel draws each chromosome's
// label centered on top of its own bar rather than beside it, so a fixed
// text color would go invisible on those bars. Standard relative-luminance
// threshold (WCAG's own formula, simplified without its gamma correction --
// overkill for a two-way black/white pick where only the ranking matters).
SYN.textColorFor = function(hex) {
    const r = parseInt(hex.slice(1, 3), 16), g = parseInt(hex.slice(3, 5), 16), b = parseInt(hex.slice(5, 7), 16);
    const luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255;
    return luminance >= 0.6 ? '#000000' : '#ffffff';
};

// reference wedges are never filtered by the min-block-size spinner (only
// ribbons are), so recoloring them only ever needs the current color count
SYN.recolorSubjectWedges = function(k, subjectSource) {
    subjectSource.data['fill_color'] = subjectSource.data['palette_index'].map((idx) => SYN.paletteColor(idx, k));
    subjectSource.change.emit();
};

// SYN.data.ribbons holds EVERY ribbon (built by SYN.buildRingLayout from
// whatever SYN.applyChainResult most recently put in SYN.data.linksByQuery/
// targetHomeologLinks/referenceHomeologLinks -- the chainer is always run at
// minBlock=3, docs/specs/hit_table.md section 4.6's least strict floor), so
// filtering by min block size, recoloring by the current color count, and
// showing/hiding self-links (is_homeolog records -- a protein hitting two
// loci within the SAME genome) all happen together here, client-side, on
// every change of any of the three controls --
// rebuilding from the same master list rather than mutating whatever
// happens to be currently displayed keeps them from stepping on each other
// (e.g. recoloring after a filter must recolor only the still-visible
// subset, not silently undo the filter or un-hide self-links).
SYN.buildOverviewRibbons = function(minScore, k, showSelfLinks, hideSynteny) {
    // self-link (is_homeolog) records are governed by showSelfLinks alone;
    // cross-genome records are governed by hideSynteny alone -- the two
    // switches compose independently rather than one overriding the other
    const records = SYN.data.ribbons.filter((r) =>
        r.score >= minScore && (r.is_homeolog ? showSelfLinks : !hideSynteny));
    return {
        xs: records.map((r) => r.xs),
        ys: records.map((r) => r.ys),
        fill_color: records.map((r) => SYN.lighten(SYN.paletteColor(r.palette_index, k))),
        alpha: records.map((r) => r.alpha),
        label: records.map((r) => r.label),
        palette_index: records.map((r) => r.palette_index),
    };
};

SYN.applyOverviewRibbons = function(minScore, k, showSelfLinks, hideSynteny, ribbonSource) {
    ribbonSource.data = SYN.buildOverviewRibbons(minScore, k, showSelfLinks, hideSynteny);
    // clear any stale ribbon-click highlight -- the rebuilt data's array
    // positions no longer correspond to whatever was selected before (a
    // filter/recolor/self-links change can add, drop, or reorder records),
    // and the fresh alphas above are already correct, so ribbonBaseAlpha
    // must be nulled BEFORE clearing .selected.indices below: that clear
    // fires SYN.applyRibbonHighlight's own listener, whose deselect path
    // only touches .data if ribbonBaseAlpha is still set -- nulling it here
    // first makes that a no-op instead of overwriting the fresh data with a
    // stale (wrong-length) snapshot from before the rebuild.
    SYN.state.ribbonBaseAlpha = null;
    ribbonSource.selected.indices = [];
    ribbonSource.change.emit();
};

// Clicking a ribbon highlights it (full opacity) and dims every other one,
// so one syntenic block can be traced by eye against the surrounding
// clutter; clicking empty ring space (or anything that clears the
// selection -- a wedge, a dotplot cell, "clear zoom") deselects and
// restores everyone's normal alpha. Rebuilds (mba/color/palette/self-links)
// also invalidate the highlight -- see SYN.applyOverviewRibbons above.
SYN.applyRibbonHighlight = function(ribbonSource) {
    const idx = ribbonSource.selected.indices;
    if (idx.length === 0) {
        if (SYN.state.ribbonBaseAlpha) {
            ribbonSource.data.alpha = SYN.state.ribbonBaseAlpha;
            SYN.state.ribbonBaseAlpha = null;
            ribbonSource.change.emit();
        }
        return;
    }
    if (!SYN.state.ribbonBaseAlpha) {
        SYN.state.ribbonBaseAlpha = ribbonSource.data.alpha.slice();
    }
    const base = SYN.state.ribbonBaseAlpha;
    const selected = new Set(idx);
    ribbonSource.data.alpha = base.map((a, i) => (selected.has(i) ? 1.0 : SYN.RIBBON_DIM_ALPHA));
    ribbonSource.change.emit();
};

// ---- ring geometry: angle math and polygon construction for the wedges
// and ribbons (point, arc, wedge, quadratic Bezier, angular chromosome
// offsets, bp-to-angle), all the way up through SYN.buildRingLayout below,
// which assembles them into a full ring for an arbitrary chromosome order --
// see module docstring for why this is the only implementation of any of
// this geometry. ----

SYN.point = function(angle, radius) {
    return [radius * Math.cos(angle), radius * Math.sin(angle)];
};

SYN.arcPoints = function(a0, a1, radius, n) {
    n = n || 48;
    if (a1 < a0) { const t = a0; a0 = a1; a1 = t; }
    const pts = [];
    for (let i = 0; i <= n; i++) { pts.push(SYN.point(a0 + (a1 - a0) * i / n, radius)); }
    return pts;
};

// Like SYN.arcPoints, but walks a0 -> a1 exactly as given, never reordered
// to a0<a1 -- SYN.arcPoints' reordering is fine for a closed wedge (outer
// arc + reversed inner arc just need consistent winding either way), but
// SYN.ribbonPolygonJS below needs to know which literal endpoint is which,
// to hand each one to the correct Bezier crossing segment.
SYN.arcPointsOrdered = function(a0, a1, radius, n) {
    n = n || 24;
    const pts = [];
    for (let i = 0; i <= n; i++) { pts.push(SYN.point(a0 + (a1 - a0) * i / n, radius)); }
    return pts;
};

SYN.wedgePolygonJS = function(a0, a1, r0, r1, n) {
    n = n || 48;
    const outer = SYN.arcPoints(a0, a1, r1, n);
    const inner = SYN.arcPoints(a0, a1, r0, n).slice().reverse();
    const poly = outer.concat(inner);
    return {xs: poly.map((p) => p[0]), ys: poly.map((p) => p[1])};
};

SYN.quadBezier = function(p0, p1, p2, n) {
    n = n || 24;
    const pts = [];
    for (let i = 0; i <= n; i++) {
        const t = i / n;
        const x = (1 - t) ** 2 * p0[0] + 2 * (1 - t) * t * p1[0] + t ** 2 * p2[0];
        const y = (1 - t) ** 2 * p0[1] + 2 * (1 - t) * t * p1[1] + t ** 2 * p2[1];
        pts.push([x, y]);
    }
    return pts;
};

SYN.ribbonPolygonJS = function(aq1, aq2, as1, as2, radius, n) {
    n = n || 24;
    // Arc across the query span -- curved along the ring's own radius, not
    // a straight chord, so it sits flush against the ideogram's own
    // curvature where the ribbon meets it (radius == innerR) -- quadratic
    // Bezier (control at the ring's center) over to the subject arc, arc
    // across the subject span, Bezier back to the start -- the bowtie
    // ribbon crossing through the middle of the ring.
    const qArc = SYN.arcPointsOrdered(aq1, aq2, radius, n);  // q1 -> q2, in order
    const sArc = SYN.arcPointsOrdered(as2, as1, radius, n);  // s2 -> s1, in order (matches traversal below)
    const origin = [0, 0];
    let poly = qArc.slice();
    poly = poly.concat(SYN.quadBezier(qArc[qArc.length - 1], origin, sArc[0], n).slice(1));
    poly = poly.concat(sArc.slice(1));
    poly = poly.concat(SYN.quadBezier(sArc[sArc.length - 1], origin, qArc[0], n).slice(1));
    return {xs: poly.map((p) => p[0]), ys: poly.map((p) => p[1])};
};

SYN.computeCircularOffsets = function(order, sizes, angleStart, angleEnd) {
    const n = order.length;
    const direction = angleEnd >= angleStart ? 1 : -1;
    const totalGap = SYN.data.chromGap * Math.max(n - 1, 0);
    const available = Math.max(Math.abs(angleEnd - angleStart) - totalGap, 0.01);
    const totalBp = order.reduce((s, name) => s + sizes[name], 0) || 1;
    const offsets = {};
    let cur = angleStart;
    for (const name of order) {
        const span = available * (sizes[name] / totalBp);
        offsets[name] = [cur, cur + direction * span];
        cur += direction * (span + SYN.data.chromGap);
    }
    return offsets;
};

SYN.bpToAngleJS = function(chrom, bpPos, sizes, offsets) {
    const [a0, a1] = offsets[chrom];
    const size = sizes[chrom];
    const frac = size ? Math.min(Math.max(bpPos / size, 0), 1) : 0;
    return a0 + frac * (a1 - a0);
};

// min/max, not a0/a1 directly: a subject/reference chromosome's angle runs
// DEcreasing with bp (see SYN.computeCircularOffsets), so a0 can be either
// side of a1 depending on which genome this gap belongs to.
SYN.widenGapAngle = function(a0, a1) {
    let lo = Math.min(a0, a1), hi = Math.max(a0, a1);
    if (hi - lo < SYN.data.minGapAngle) {
        const mid = (lo + hi) / 2;
        lo = mid - SYN.data.minGapAngle / 2;
        hi = mid + SYN.data.minGapAngle / 2;
    }
    return [lo, hi];
};

// Zoom-panel gap ticks -- a purely visual visibility fix, same idea as
// SYN.widenGapAngle above
// (a gap is typically a few hundred bp against a multi-Mb bar, so drawn at
// true width most would render sub-pixel) but in linear bp/px space instead
// of angular: refSpan is the caller's own best estimate of the panel's full
// data-space width (pivotSize where that's all that's known yet, xMax once
// it's computed -- see each call site), converted through the panel's fixed
// SYN.data.detailFigWidth px size and the same x_range padding build_page()
// uses (xMax * 1.08) to get an approximate bp-per-px, then floored to
// SYN.GAP_MIN_PX pixels' worth of bp. Small on purpose: a gap tick spans the
// bar's FULL height (no vertical inset -- see each addGapTicks call site's
// own gapYs), so it should read as a thin needle across that height, not a
// short wide block -- only the width needs flooring, not the height.
SYN.GAP_MIN_PX = 1.5;
SYN.widenGapBp = function(start, end, refSpan) {
    const minBp = (refSpan * 1.08 / SYN.data.detailFigWidth) * SYN.GAP_MIN_PX;
    if (end - start < minBp) {
        const mid = (start + end) / 2;
        return [mid - minBp / 2, mid + minBp / 2];
    }
    return [start, end];
};

// Rebuilds every ring-native source (query/subject wedges, chromosome
// labels, and the ribbons master list SYN.data.ribbons itself -- NOT just
// whatever subset is currently filtered into r_src, so a later mba/color/
// self-links change still has correct fresh geometry to filter from) for an
// arbitrary (queryOrder, subjectOrder). The only place any of this geometry
// is built -- SYN.init (see build_page()) calls this for the page's very
// first render too, through the same path as every later reorder.
// Ribbon master list (SYN.data.ribbons) for the given ring offsets -- split out
// of SYN.buildRingLayout so a re-chain (SYN.applyChainResult) can rebuild just
// the ribbons: the wedges, labels and assembly-gap marks depend only on the
// chromosome order, never on the chain result, and on a gappy draft
// assembly the gap marks alone are the costliest part of a full layout.
SYN.buildRibbonRecords = function(queryOffsets, subjectOffsets) {
    // Each of the three groups below is sorted by score ascending (stable)
    // before its ribbons are pushed, so within each group a higher-scoring
    // ribbon always ends up later in the array -- and therefore drawn on
    // top, Bokeh's own patches() convention -- regardless of what order
    // SYN.data.linksByQuery/targetHomeologLinks/referenceHomeologLinks
    // themselves happen to iterate in, and regardless of how many times the
    // ring has been reordered since page load.
    const ribbons = [];
    const allCross = [].concat(...Object.values(SYN.data.linksByQuery)).sort((a, b) => a.score - b.score);
    for (const l of allCross) {
        // either endpoint's chromosome may have been dropped by the min-
        // length filter (see SYN.filterBySize) -- bpToAngleJS would throw
        // destructuring offsets[chrom] === undefined, so skip rather than
        // draw a ribbon to a chromosome no longer on the ring at all
        if (!(l.q_chrom in queryOffsets) || !(l.s_chrom in subjectOffsets)) { continue; }
        const aq1 = SYN.bpToAngleJS(l.q_chrom, l.q_start, SYN.data.querySizes, queryOffsets);
        const aq2 = SYN.bpToAngleJS(l.q_chrom, l.q_end, SYN.data.querySizes, queryOffsets);
        const as1 = SYN.bpToAngleJS(l.s_chrom, l.s_start, SYN.data.subjectSizes, subjectOffsets);
        const as2 = SYN.bpToAngleJS(l.s_chrom, l.s_end, SYN.data.subjectSizes, subjectOffsets);
        const poly = SYN.ribbonPolygonJS(aq1, aq2, as1, as2, SYN.data.linkR);
        ribbons.push({
            xs: poly.xs, ys: poly.ys, palette_index: SYN.data.subjectColorIndex[l.s_chrom],
            alpha: 0.25 + 0.55 * (l.score / SYN.data.maxScore), label: SYN.formatLinkLabel(l),
            score: l.score, is_homeolog: false,
        });
    }
    const targetHomeologLinks = SYN.data.targetHomeologLinks.slice().sort((a, b) => a.score - b.score);
    for (const l of targetHomeologLinks) {
        if (!(l.q_chrom in queryOffsets) || !(l.s_chrom in queryOffsets)) { continue; }
        const aq1 = SYN.bpToAngleJS(l.q_chrom, l.q_start, SYN.data.querySizes, queryOffsets);
        const aq2 = SYN.bpToAngleJS(l.q_chrom, l.q_end, SYN.data.querySizes, queryOffsets);
        const as1 = SYN.bpToAngleJS(l.s_chrom, l.s_start, SYN.data.querySizes, queryOffsets);
        const as2 = SYN.bpToAngleJS(l.s_chrom, l.s_end, SYN.data.querySizes, queryOffsets);
        const poly = SYN.ribbonPolygonJS(aq1, aq2, as1, as2, SYN.data.linkR);
        ribbons.push({
            xs: poly.xs, ys: poly.ys, palette_index: SYN.data.queryColorIndex[l.q_chrom],
            alpha: 0.25 + 0.55 * (l.score / SYN.data.maxScore),
            label: SYN.formatLinkLabel(l, 'Target', 'Target') + ' (homeolog)',
            score: l.score, is_homeolog: true,
        });
    }
    const referenceHomeologLinks = SYN.data.referenceHomeologLinks.slice().sort((a, b) => a.score - b.score);
    for (const l of referenceHomeologLinks) {
        if (!(l.q_chrom in subjectOffsets) || !(l.s_chrom in subjectOffsets)) { continue; }
        const aq1 = SYN.bpToAngleJS(l.q_chrom, l.q_start, SYN.data.subjectSizes, subjectOffsets);
        const aq2 = SYN.bpToAngleJS(l.q_chrom, l.q_end, SYN.data.subjectSizes, subjectOffsets);
        const as1 = SYN.bpToAngleJS(l.s_chrom, l.s_start, SYN.data.subjectSizes, subjectOffsets);
        const as2 = SYN.bpToAngleJS(l.s_chrom, l.s_end, SYN.data.subjectSizes, subjectOffsets);
        const poly = SYN.ribbonPolygonJS(aq1, aq2, as1, as2, SYN.data.linkR);
        ribbons.push({
            xs: poly.xs, ys: poly.ys, palette_index: SYN.data.subjectColorIndex[l.q_chrom],
            alpha: 0.25 + 0.55 * (l.score / SYN.data.maxScore),
            label: SYN.formatLinkLabel(l, 'Reference', 'Reference') + ' (homeolog)',
            score: l.score, is_homeolog: true,
        });
    }
    return ribbons;
};

SYN.buildRingLayout = function(queryOrder, subjectOrder, k) {
    const groupGap = SYN.data.groupGap;
    const subjectOffsets = SYN.computeCircularOffsets(
        subjectOrder, SYN.data.subjectSizes, Math.PI - groupGap, groupGap);
    const queryOffsets = SYN.computeCircularOffsets(
        queryOrder, SYN.data.querySizes, Math.PI + groupGap, 2 * Math.PI - groupGap);

    const q = {xs: [], ys: [], fill_color: [], name: [], size_label: [], group: []};
    for (const name of queryOrder) {
        const [a0, a1] = queryOffsets[name];
        const poly = SYN.wedgePolygonJS(a0, a1, SYN.data.innerR, SYN.data.outerR);
        q.xs.push(poly.xs); q.ys.push(poly.ys);
        q.fill_color.push(SYN.data.targetGrey); q.name.push(name);
        q.size_label.push(`${(SYN.data.querySizes[name] / 1e6).toFixed(2)} Mb`);
        q.group.push('target');
    }

    const s = {xs: [], ys: [], fill_color: [], name: [], size_label: [], group: [], palette_index: []};
    for (const name of subjectOrder) {
        const [a0, a1] = subjectOffsets[name];
        const poly = SYN.wedgePolygonJS(a0, a1, SYN.data.innerR, SYN.data.outerR);
        const pIdx = SYN.data.subjectColorIndex[name];
        s.xs.push(poly.xs); s.ys.push(poly.ys);
        s.fill_color.push(SYN.paletteColor(pIdx, k));
        s.name.push(name);
        s.size_label.push(`${(SYN.data.subjectSizes[name] / 1e6).toFixed(2)} Mb`);
        s.group.push('reference'); s.palette_index.push(pIdx);
    }

    const labelR = SYN.data.outerR + 0.14;
    const label = {x: [], y: [], text: []};
    for (const name of queryOrder) {
        const [a0, a1] = queryOffsets[name]; const mid = (a0 + a1) / 2;
        label.x.push(labelR * Math.cos(mid)); label.y.push(labelR * Math.sin(mid)); label.text.push(name);
    }
    for (const name of subjectOrder) {
        const [a0, a1] = subjectOffsets[name]; const mid = (a0 + a1) / 2;
        label.x.push(labelR * Math.cos(mid)); label.y.push(labelR * Math.sin(mid)); label.text.push(name);
    }

    const ribbons = SYN.buildRibbonRecords(queryOffsets, subjectOffsets);

    // Assembly-gap wedges for the current order -- always computed here
    // regardless of show_gaps_toggle's current state (mirrors ribbons above:
    // the toggle only ever filters what SYN.applyGapVisibility copies OUT of
    // this master list into gap_src, same pattern as
    // SYN.buildOverviewRibbons/SYN.applyOverviewRibbons for self-links).
    const gaps = SYN.buildGapRecords(queryOffsets, subjectOffsets, SYN.GAP_ARC_SEGMENTS);

    return {q, s, label, ribbons, gaps, queryOffsets, subjectOffsets};
};

// A gap wedge is at most a fraction of a degree wide (see MIN_GAP_ANGLE's
// own Python comment), so SYN.wedgePolygonJS's 48-segment default arc is
// wasted precision -- 2 segments (a straight line across that width) reads
// identically to a true arc at this scale, and cuts every gap wedge from
// 98 vertices to 8. A parameter, not a hardcoded literal inside
// buildGapRecords, so the page_geometry.mjs test can call it with 48 to
// reproduce the old page's full-resolution gap geometry for comparison.
SYN.GAP_ARC_SEGMENTS = 2;

// Ring-wedge geometry for a genome's own assembly gaps -- same annular-
// sector shape as a chromosome wedge, same full innerR-outerR radial span
// too (SYN.wedgePolygonJS's own r0/r1 args below, no radial inset -- a gap
// marker is exactly as tall as the ideogram it marks), just angularly
// widened to SYN.data.minGapAngle if narrower (see that constant's own
// Python comment) so it reads as a thin needle across that height rather
// than a sub-pixel sliver. Pulled out of SYN.buildRingLayout so
// page_geometry.mjs can call it directly at nSegments=48 to compare against
// the old page's baked (also 48-segment) gap wedges.
SYN.buildGapRecords = function(queryOffsets, subjectOffsets, nSegments) {
    const gaps = [];
    for (const chrom in SYN.data.targetGapsByChrom) {
        if (!(chrom in queryOffsets)) { continue; }
        for (const g of SYN.data.targetGapsByChrom[chrom]) {
            const a0 = SYN.bpToAngleJS(chrom, g.start, SYN.data.querySizes, queryOffsets);
            const a1 = SYN.bpToAngleJS(chrom, g.end, SYN.data.querySizes, queryOffsets);
            const [lo, hi] = SYN.widenGapAngle(a0, a1);
            const poly = SYN.wedgePolygonJS(lo, hi, SYN.data.innerR, SYN.data.outerR, nSegments);
            gaps.push({xs: poly.xs, ys: poly.ys, label: SYN.formatGapLabel(chrom, g.start, g.end)});
        }
    }
    for (const chrom in SYN.data.referenceGapsByChrom) {
        if (!(chrom in subjectOffsets)) { continue; }
        for (const g of SYN.data.referenceGapsByChrom[chrom]) {
            const a0 = SYN.bpToAngleJS(chrom, g.start, SYN.data.subjectSizes, subjectOffsets);
            const a1 = SYN.bpToAngleJS(chrom, g.end, SYN.data.subjectSizes, subjectOffsets);
            const [lo, hi] = SYN.widenGapAngle(a0, a1);
            const poly = SYN.wedgePolygonJS(lo, hi, SYN.data.innerR, SYN.data.outerR, nSegments);
            gaps.push({xs: poly.xs, ys: poly.ys, label: SYN.formatGapLabel(chrom, g.start, g.end)});
        }
    }
    return gaps;
};

// Applies a full ring reorder: rebuilds wedge/label sources directly (never
// filtered by anything), replaces the SYN.data.ribbons MASTER list with
// fresh geometry (so mba/color/self-links/hide-synteny keep working
// correctly against the new order afterward), then re-runs the current
// filter state on top of it via SYN.applyOverviewRibbons. Clears wedge/
// ribbon selections -- stale indices from before the rebuild would point at
// the wrong chromosome/ribbon after the arrays are replaced (same reasoning
// as SYN.applyOverviewRibbons's own selection-clearing).
SYN.applyRingLayout = function(queryOrder, subjectOrder, sources) {
    const layout = SYN.buildRingLayout(queryOrder, subjectOrder, sources.k);
    sources.querySource.data = layout.q;
    sources.subjectSource.data = layout.s;
    sources.labelSource.data = layout.label;
    SYN.data.ribbons = layout.ribbons;
    SYN.data.gapRecords = layout.gaps;
    SYN.state.ringQueryOffsets = layout.queryOffsets;
    SYN.state.ringSubjectOffsets = layout.subjectOffsets;
    sources.querySource.selected.indices = [];
    sources.subjectSource.selected.indices = [];
    sources.querySource.change.emit();
    sources.subjectSource.change.emit();
    sources.labelSource.change.emit();
    SYN.applyOverviewRibbons(sources.minScore, sources.k, sources.showSelfLinks,
                              sources.hideSynteny, sources.ribbonSource);
    SYN.applyGapVisibility(sources.gapSource);
};

// SYN.data.gapRecords holds the CURRENT ring order's full wedge-polygon
// list (every gap, replaced wholesale on a reorder -- see SYN.applyRingLayout
// above), so this only ever needs to copy it into gap_src, or clear it, per
// SYN.state.showGaps -- same "rebuild from a master list" idiom as
// SYN.applyOverviewRibbons, just with no score/self-links filtering to do.
SYN.applyGapVisibility = function(gapSource) {
    const records = SYN.state.showGaps ? SYN.data.gapRecords : [];
    gapSource.data = {
        xs: records.map((r) => r.xs), ys: records.map((r) => r.ys), label: records.map((r) => r.label),
    };
    gapSource.change.emit();
};

// The dotplot's chromosome order is not fixed the way the ring's is -- the
// "order by similarity" toggle (see build_page()) can replace it wholesale,
// so every dotplot source that depends on chromosome order (rulers, grid,
// labels, cells, and the block segments themselves) is rebuilt from raw
// per-link data on every redraw rather than only once at page load, unlike
// SYN.buildOverviewRibbons above (the ring's order never changes, so its
// precomputed SYN.data.ribbons records stay valid for the page's whole
// lifetime and only need filtering/recoloring, never rebuilding from
// scratch). SYN.state.dpQueryOffsets/dpSubjectOffsets hold whichever order
// is currently active; this and buildDotplotLayout below are the only two
// places that read them.
SYN.buildDotplotSegmentsForLayout = function(queryOffsets, subjectOffsets, minScore, k) {
    const xs = [], ys = [], lineColor = [], alpha = [], label = [], paletteIndex = [];
    const maxScore = SYN.data.maxScore;
    for (const qName in SYN.data.linksByQuery) {
        for (const l of SYN.data.linksByQuery[qName]) {
            if (l.score < minScore) { continue; }
            // either chromosome may be below the min-length filter (see
            // SYN.filterBySize) and absent from this order's offsets
            if (!(l.q_chrom in queryOffsets) || !(l.s_chrom in subjectOffsets)) { continue; }
            const x0 = queryOffsets[l.q_chrom] + l.q_start;
            const x1 = queryOffsets[l.q_chrom] + l.q_end;
            const yLo = subjectOffsets[l.s_chrom] + l.s_start;
            const yHi = subjectOffsets[l.s_chrom] + l.s_end;
            // '-' (inverted) blocks run the reference span the opposite
            // direction from the target span -- the same convention
            // MCScanX/D-GENIES-style dotplots use to make inversions
            // visually distinct from collinear blocks
            const [y0, y1] = l.orientation !== '-' ? [yLo, yHi] : [yHi, yLo];
            xs.push([x0, x1]); ys.push([y0, y1]);
            const pIdx = SYN.data.subjectColorIndex[l.s_chrom];
            lineColor.push(SYN.lighten(SYN.paletteColor(pIdx, k)));
            alpha.push(0.25 + 0.55 * (l.score / maxScore));
            label.push(SYN.formatLinkLabel(l));
            paletteIndex.push(pIdx);
        }
    }
    return {xs, ys, line_color: lineColor, alpha, label, palette_index: paletteIndex};
};

SYN.applyDotplotSegmentsForCurrentLayout = function(minScore, k, segmentSource) {
    segmentSource.data = SYN.buildDotplotSegmentsForLayout(SYN.state.dpQueryOffsets, SYN.state.dpSubjectOffsets, minScore, k);
    segmentSource.change.emit();
};

SYN.computeOffsets = function(order, sizes) {
    const offsets = {};
    let cur = 0;
    for (const name of order) { offsets[name] = cur; cur += sizes[name]; }
    return offsets;
};

// The non-similarity order: target reads left-to-right, reference
// bottom-to-top, both in whichever of {size, natural (FASTA)} order
// SYN.state.orderBySize currently selects (see size_order_toggle) -- no
// reversal (an earlier version reversed the reference list, putting the
// first name at the top -- flipped per explicit request): SYN.computeOffsets
// gives the first name in its list the lowest -- bottommost -- position,
// same convention the x-axis already uses for the target's own offsets.
// Named "natural" for what it's the alternative TO (similarity order, see
// SYN.computeSimilarityOrder below) -- not because it always resolves to
// FASTA order anymore.
SYN.dpNaturalOrder = function() {
    const bySize = SYN.state.orderBySize;
    return {
        queryOrder: (bySize ? SYN.data.queryNames : SYN.data.queryNamesNatural).slice(),
        subjectOrder: (bySize ? SYN.data.subjectNames : SYN.data.subjectNamesNatural).slice(),
    };
};

// Drops any chromosome shorter than SYN.state.minSeqSize (see
// min_seq_size_spinner in build_page()) from an order list -- composes with
// whichever ordering is currently active rather than being its own separate
// order, since min_seq_size_spinner's own callback (and every other control
// that can change the active order) always re-derives its order through
// SYN.ringOrderFor/SYN.dpOrderFor below, never the raw name lists directly.
SYN.filterBySize = function(order, sizes) {
    return order.filter((name) => sizes[name] >= SYN.state.minSeqSize);
};

// The ring's currently-active order (size vs. natural -- it has no
// similarity concept, see SYN.buildRingLayout's own comment), length-filtered.
SYN.ringOrderFor = function(bySize) {
    return {
        queryOrder: SYN.filterBySize(bySize ? SYN.data.queryNames : SYN.data.queryNamesNatural, SYN.data.querySizes),
        subjectOrder: SYN.filterBySize(bySize ? SYN.data.subjectNames : SYN.data.subjectNamesNatural, SYN.data.subjectSizes),
    };
};

// The dotplot's currently-active order (similarity vs. SYN.dpNaturalOrder's
// size/natural pair, mirroring order_toggle's own ternary), length-filtered.
SYN.dpOrderFor = function(useSimilarity) {
    const order = useSimilarity ? SYN.computeSimilarityOrder() : SYN.dpNaturalOrder();
    return {
        queryOrder: SYN.filterBySize(order.queryOrder, SYN.data.querySizes),
        subjectOrder: SYN.filterBySize(order.subjectOrder, SYN.data.subjectSizes),
    };
};

// Reorders both axes to make shared synteny read as a diagonal: each
// chromosome is placed near the average position of whatever it shares
// anchors with on the other axis (weighted by score), alternating which
// axis is being reordered against the other's current order a few times so
// the two converge together rather than one being optimized against a
// now-stale copy of the other. A chromosome with no cross-genome anchors at
// all (score-weighted total of 0) sorts to the end rather than landing at
// an arbitrary position -- there's no similarity signal to place it by.
SYN.computeSimilarityOrder = function() {
    const natural = SYN.dpNaturalOrder();
    let queryOrder = natural.queryOrder, subjectOrder = natural.subjectOrder;

    const weight = (q, s) => {
        let total = 0;
        for (const l of (SYN.data.linksByQuery[q] || [])) { if (l.s_chrom === s) { total += l.score; } }
        return total;
    };

    const reorder = (names, otherOrder, weightOf) => {
        const rankOf = {};
        otherOrder.forEach((name, i) => { rankOf[name] = i; });
        return names
            .map((name, i) => {
                let totalW = 0, weightedRank = 0;
                for (const other of otherOrder) {
                    const w = weightOf(name, other);
                    if (w > 0) { totalW += w; weightedRank += w * rankOf[other]; }
                }
                return {name, i, avgRank: totalW > 0 ? weightedRank / totalW : Infinity};
            })
            // the index tiebreak keeps a stable, deterministic order among
            // chromosomes with no signal, instead of leaving their relative
            // order up to the sort algorithm's whim
            .sort((a, b) => (a.avgRank - b.avgRank) || (a.i - b.i))
            .map((r) => r.name);
    };

    for (let iter = 0; iter < 3; iter++) {
        queryOrder = reorder(queryOrder, subjectOrder, weight);
        subjectOrder = reorder(subjectOrder, queryOrder, (s, q) => weight(q, s));
    }
    return {queryOrder, subjectOrder};
};

// Full rebuild of every order-dependent dotplot source for an arbitrary
// (queryOrder, subjectOrder) pair -- the only implementation of any of this
// geometry (see module docstring), called for every reorder, recolor, or
// refilter, and for the page's very first render too (see SYN.init).
// rx/ry/totalX/totalY are recomputed on every call rather than assumed
// constant: a pure reorder never changes them, but the min-length filter
// (see SYN.filterBySize) can DROP chromosomes between calls, which does
// change the total span.
SYN.buildDotplotLayout = function(queryOrder, subjectOrder, k) {
    // the `|| 1` floor mirrors SYN.computeCircularOffsets' own totalBp
    // fallback -- keeps the figure/ruler math sane (not NaN/zero-width) in
    // the degenerate case where the filter has emptied one whole side
    const totalX = queryOrder.reduce((sum, name) => sum + SYN.data.querySizes[name], 0) || 1;
    const totalY = subjectOrder.reduce((sum, name) => sum + SYN.data.subjectSizes[name], 0) || 1;
    const rx = totalX * SYN.data.dpRulerFrac, ry = totalY * SYN.data.dpRulerFrac;
    const queryOffsets = SYN.computeOffsets(queryOrder, SYN.data.querySizes);
    const subjectOffsets = SYN.computeOffsets(subjectOrder, SYN.data.subjectSizes);

    const qXs = [], qYs = [], qName = [], qSize = [], qGroup = [];
    for (const name of queryOrder) {
        const size = SYN.data.querySizes[name], x0 = queryOffsets[name];
        qXs.push([x0, x0 + size, x0 + size, x0]); qYs.push([-ry, -ry, 0, 0]);
        qName.push(name); qSize.push((size / 1e6).toFixed(2) + ' Mb'); qGroup.push('target');
    }
    const sXs = [], sYs = [], sName = [], sSize = [], sGroup = [], sIdx = [];
    for (const name of subjectOrder) {
        const size = SYN.data.subjectSizes[name], y0 = subjectOffsets[name];
        sXs.push([-rx, -rx, 0, 0]); sYs.push([y0, y0 + size, y0 + size, y0]);
        sName.push(name); sSize.push((size / 1e6).toFixed(2) + ' Mb'); sGroup.push('reference');
        sIdx.push(SYN.data.subjectColorIndex[name]);
    }

    const gridXs = [], gridYs = [];
    for (const name of queryOrder) {
        const x0 = queryOffsets[name];
        gridXs.push([x0, x0]); gridYs.push([-ry, totalY]);
    }
    // the min-length filter (see SYN.filterBySize) can empty one whole side
    // -- nothing to draw a trailing boundary line for in that case
    if (queryOrder.length) {
        const lastQ = queryOrder[queryOrder.length - 1];
        gridXs.push([queryOffsets[lastQ] + SYN.data.querySizes[lastQ], queryOffsets[lastQ] + SYN.data.querySizes[lastQ]]);
        gridYs.push([-ry, totalY]);
    }
    for (const name of subjectOrder) {
        const y0 = subjectOffsets[name];
        gridXs.push([-rx, totalX]); gridYs.push([y0, y0]);
    }
    if (subjectOrder.length) {
        const lastS = subjectOrder[subjectOrder.length - 1];
        gridXs.push([-rx, totalX]);
        gridYs.push([subjectOffsets[lastS] + SYN.data.subjectSizes[lastS], subjectOffsets[lastS] + SYN.data.subjectSizes[lastS]]);
    }

    const qLabelX = [], qLabelY = [], qLabelText = [];
    for (const name of queryOrder) {
        const size = SYN.data.querySizes[name], x0 = queryOffsets[name];
        qLabelX.push(x0 + size / 2); qLabelY.push(-ry * 1.7); qLabelText.push(name);
    }
    const sLabelX = [], sLabelY = [], sLabelText = [];
    for (const name of subjectOrder) {
        const size = SYN.data.subjectSizes[name], y0 = subjectOffsets[name];
        sLabelX.push(-rx * 1.7); sLabelY.push(y0 + size / 2); sLabelText.push(name);
    }

    return {
        queryOffsets, subjectOffsets, totalX, totalY, rx, ry,
        queryRuler: {xs: qXs, ys: qYs, fill_color: queryOrder.map(() => SYN.data.targetGrey),
                     name: qName, size_label: qSize, group: qGroup},
        subjectRuler: {xs: sXs, ys: sYs, fill_color: sIdx.map((idx) => SYN.paletteColor(idx, k)),
                       name: sName, size_label: sSize, group: sGroup, palette_index: sIdx},
        grid: {xs: gridXs, ys: gridYs},
        queryLabels: {x: qLabelX, y: qLabelY, text: qLabelText},
        subjectLabels: {x: sLabelX, y: sLabelY, text: sLabelText},
    };
};

// Applies a full reorder: rebuilds every order-dependent dotplot source (via
// buildDotplotLayout) and updates SYN.state's offsets so every later
// recolor/refilter (color spinner, palette dropdown, min-block-size) draws
// against the new order without needing to know a reorder happened.
// sources.fig, if given, is the dotplot figure itself -- its x_range/y_range
// are refreshed to the new totalX/totalY so a min-length filter (see
// SYN.filterBySize) that shrinks the visible chromosome set doesn't leave a
// stale range showing dead space beyond wherever the plot now actually ends;
// a no-op for a pure reorder, since the total is unchanged then.
SYN.applyDotplotOrder = function(queryOrder, subjectOrder, k, minScore, sources) {
    const layout = SYN.buildDotplotLayout(queryOrder, subjectOrder, k);
    SYN.state.dpQueryOffsets = layout.queryOffsets;
    SYN.state.dpSubjectOffsets = layout.subjectOffsets;
    // the order arrays themselves (not just their offsets) -- SYN.dpCellAt
    // binary-searches these to hit-test a dotplot tap without a per-cell
    // renderer (see module docstring on why that layer is gone)
    SYN.state.dpQueryOrder = queryOrder;
    SYN.state.dpSubjectOrder = subjectOrder;
    SYN.state.dpTotalX = layout.totalX; SYN.state.dpTotalY = layout.totalY;
    SYN.state.dpRx = layout.rx; SYN.state.dpRy = layout.ry;
    sources.query.data = layout.queryRuler; sources.query.change.emit();
    sources.subject.data = layout.subjectRuler; sources.subject.change.emit();
    sources.grid.data = layout.grid; sources.grid.change.emit();
    sources.queryLabel.data = layout.queryLabels; sources.queryLabel.change.emit();
    sources.subjectLabel.data = layout.subjectLabels; sources.subjectLabel.change.emit();
    if (sources.fig) {
        // *3 (not just enough room for the ruler strip itself) leaves space
        // for the chromosome-name labels drawn just outside each ruler (see
        // buildDotplotLayout's own qLabelX/sLabelX -- text_align='right' for
        // the reference labels means they extend further left from their
        // anchor, so the range has to clear that too, not just the ruler).
        // Confirmed against real data: *2.2 and *2.6 clip a 5-character name
        // like "chr13" down to "13"/"hr13"; *3 leaves every label fully
        // visible with a little room to spare while still keeping the dead
        // margin well under half of what a much larger multiplier would waste.
        sources.fig.x_range.start = -layout.rx * 3; sources.fig.x_range.end = layout.totalX * 1.02;
        sources.fig.y_range.start = -layout.ry * 3; sources.fig.y_range.end = layout.totalY * 1.02;
    }
    SYN.applyDotplotSegmentsForCurrentLayout(minScore, k, sources.segment);
    if (sources.gap) { SYN.applyDotplotGapVisibility(sources.gap); }
};

// Binary search for the one name in `order` (already position-sorted along
// this axis by SYN.computeOffsets) whose [offset, offset+size) interval
// contains pos -- the dotplot lays chromosomes flush end to end with no gap
// (see SYN.computeOffsets), so this always resolves to exactly one name
// inside [0, total), never a gap between two names. Returns null for a pos
// outside every interval (caller's job to bounds-check against the total
// first -- see SYN.dpCellAt).
SYN.chromAtOffset = function(order, offsets, sizes, pos) {
    let lo = 0, hi = order.length - 1;
    while (lo <= hi) {
        const mid = (lo + hi) >> 1;
        const name = order[mid];
        const start = offsets[name], end = start + sizes[name];
        if (pos < start) { hi = mid - 1; }
        else if (pos >= end) { lo = mid + 1; }
        else { return name; }
    }
    return null;
};

// Replaces the old per-(target, reference) invisible grid-cell renderer
// (48,841 hidden polygons for axolotl -- see module docstring) with a plain
// hit test: two binary searches, one per axis, against whichever dotplot
// order is currently active (SYN.state.dpQueryOrder/dpSubjectOrder, kept
// current by SYN.applyDotplotOrder). Returns null outside the grid entirely
// (the ruler strips at negative x/y, and the dead margin beyond either axis'
// total) -- see dotplot_fig's own Tap handler in build_page(), the only
// caller, for why that's exactly the cases where a plain tap should fall
// through to the rulers' own TapTool-driven selection instead.
SYN.dpCellAt = function(x, y) {
    if (x < 0 || y < 0 || x >= SYN.state.dpTotalX || y >= SYN.state.dpTotalY) { return null; }
    const targetName = SYN.chromAtOffset(SYN.state.dpQueryOrder, SYN.state.dpQueryOffsets, SYN.data.querySizes, x);
    const subjectName = SYN.chromAtOffset(SYN.state.dpSubjectOrder, SYN.state.dpSubjectOffsets, SYN.data.subjectSizes, y);
    if (targetName === null || subjectName === null) { return null; }
    return {targetName, subjectName};
};

// Dotplot companion to SYN.buildRingLayout's gap wedges (nice-to-have, same
// "Show gaps" switch) -- thin lines instead of a polygon per gap, same
// bounds convention as SYN.buildDotplotLayout's own boundary gridlines --
// entirely client-side, reading SYN.state.dpQueryOffsets/
// dpSubjectOffsets/dpRx/dpRy/dpTotalX/dpTotalY directly (always current --
// SYN.applyDotplotOrder sets all of them before this ever runs) rather than
// taking offsets as parameters, so show_gaps_toggle's own callback can call
// this without needing to know which order (or min-length filter, see
// SYN.filterBySize) is currently active.
SYN.buildDotplotGapLines = function() {
    const rx = SYN.state.dpRx, ry = SYN.state.dpRy;
    const totalX = SYN.state.dpTotalX, totalY = SYN.state.dpTotalY;
    const xs = [], ys = [], label = [];
    for (const chrom in SYN.data.targetGapsByChrom) {
        // dropped by the min-length filter -- no offset to place this at
        if (!(chrom in SYN.state.dpQueryOffsets)) { continue; }
        for (const g of SYN.data.targetGapsByChrom[chrom]) {
            const mid = SYN.state.dpQueryOffsets[chrom] + (g.start + g.end) / 2;
            xs.push([mid, mid]); ys.push([-ry, totalY]);
            label.push(SYN.formatGapLabel(chrom, g.start, g.end));
        }
    }
    for (const chrom in SYN.data.referenceGapsByChrom) {
        if (!(chrom in SYN.state.dpSubjectOffsets)) { continue; }
        for (const g of SYN.data.referenceGapsByChrom[chrom]) {
            const mid = SYN.state.dpSubjectOffsets[chrom] + (g.start + g.end) / 2;
            xs.push([-rx, totalX]); ys.push([mid, mid]);
            label.push(SYN.formatGapLabel(chrom, g.start, g.end));
        }
    }
    return {xs, ys, label};
};

SYN.applyDotplotGapVisibility = function(gapSource) {
    gapSource.data = SYN.state.showGaps ? SYN.buildDotplotGapLines() : {xs: [], ys: [], label: []};
    gapSource.change.emit();
};

SYN.buildOffsetsLinear = function(chroms, gap) {
    const offsets = {};
    let cum = 0;
    for (const [name, size] of chroms) {
        offsets[name] = cum;
        cum += size + gap;
    }
    return [offsets, cum - gap];
};

// 'target' -> SYN.data.targetGapsByChrom, 'reference' -> ...GapsByChrom --
// both the zoom panel builders below need this same per-chromosome lookup
// (as opposed to the ring's flat wedge-polygon list, see SYN.buildRingLayout)
// to know which genome's gap list a given bar's chromosome name belongs to,
// since target and reference chromosomes can share names (see
// Dataset.__init__'s own comment on this).
SYN.gapsForChrom = function(role, name) {
    return (role === 'target' ? SYN.data.targetGapsByChrom : SYN.data.referenceGapsByChrom)[name] || [];
};

SYN.emptyDetail = function(title) {
    return {
        bars: {xs: [], ys: [], fill_color: []},
        ribbons: {xs: [], ys: [], fill_color: [], alpha: [], label: []},
        labels: {x: [], y: [], text: [], color: []},
        gaps: {xs: [], ys: [], label: []},
        x_range: [-1, 1], y_range: [-0.5, 1.5],
        title: title || 'Click any chromosome wedge on the left to zoom in.',
    };
};

// pivotSide is 'subject' (reference) or 'query' (target) -- whichever side
// was clicked. That chromosome is drawn alone on its own row; the other
// side's chromosomes involved in its links are packed on the opposite row.
// Reference-associated bars/ribbons always end up on the top row, target-
// associated ones on the bottom row, regardless of which side is the pivot.
SYN.buildDetailData = function(pivotSide, pivotName, k, minScore) {
    if (!pivotName) { return SYN.emptyDetail(); }
    const isPivotSubject = pivotSide === 'subject';
    const links = ((isPivotSubject ? SYN.data.linksByReference[pivotName] : SYN.data.linksByQuery[pivotName]) || [])
        .filter((l) => l.score >= minScore);
    const pivotSize = isPivotSubject ? SYN.data.subjectSizes[pivotName] : SYN.data.querySizes[pivotName];
    const pivotColor = isPivotSubject ? SYN.paletteColor(SYN.data.subjectColorIndex[pivotName], k) : SYN.data.targetGrey;
    const pivotKind = isPivotSubject ? SYN.state.referenceLabel : SYN.state.targetLabel;
    const otherKind = isPivotSubject ? SYN.state.targetLabel : SYN.state.referenceLabel;

    const barH = 0.32, topY = 1.2, botY = 0.0;
    const pivotRowY = isPivotSubject ? topY : botY;
    const otherRowY = isPivotSubject ? botY : topY;
    const pivotRole = isPivotSubject ? 'reference' : 'target';
    const otherRole = isPivotSubject ? 'target' : 'reference';

    // gap ticks for one bar -- same thin-rect-over-the-bar idea as the
    // ring's wedge markers, just linear: x0 is the bar's own left edge (0
    // for the pivot, which is never offset; otherOffsets[name] for a packed
    // "other" bar), rowY/barH match whichever row that bar is drawn on
    // (full bar height, no vertical inset -- a gap tick is exactly as tall
    // as the bar it marks, same as the ring's gap wedges), and refSpan is
    // the caller's best-known estimate of the panel's overall bp span (see
    // SYN.widenGapBp).
    const gapXs = [], gapYs = [], gapLabel = [];
    const addGapTicks = (role, name, x0, rowY, refSpan) => {
        for (const g of SYN.gapsForChrom(role, name)) {
            const [gs, ge] = SYN.widenGapBp(g.start, g.end, refSpan);
            gapXs.push([x0 + gs, x0 + ge, x0 + ge, x0 + gs]);
            gapYs.push([rowY, rowY, rowY + barH, rowY + barH]);
            gapLabel.push(SYN.formatGapLabel(name, g.start, g.end));
        }
    };

    if (links.length === 0) {
        const d = SYN.emptyDetail(`${pivotName}: no links found`);
        d.bars = {xs: [[0, pivotSize, pivotSize, 0]],
                  ys: [[pivotRowY, pivotRowY, pivotRowY + barH, pivotRowY + barH]],
                  fill_color: [pivotColor]};
        d.labels = {x: [pivotSize / 2], y: [pivotRowY + barH / 2], text: [pivotName],
                    color: [SYN.textColorFor(pivotColor)]};
        addGapTicks(pivotRole, pivotName, 0, pivotRowY, pivotSize);
        d.gaps = {xs: gapXs, ys: gapYs, label: gapLabel};
        d.x_range = [-pivotSize * 0.03, pivotSize * 1.05];
        return d;
    }

    const otherName = (l) => isPivotSubject ? l.q_chrom : l.s_chrom;
    const otherStart = (l) => isPivotSubject ? l.q_start : l.s_start;
    const otherEnd = (l) => isPivotSubject ? l.q_end : l.s_end;
    const pivotStart = (l) => isPivotSubject ? l.s_start : l.q_start;
    const pivotEnd = (l) => isPivotSubject ? l.s_end : l.q_end;
    const otherSizes = isPivotSubject ? SYN.data.querySizes : SYN.data.subjectSizes;

    const otherNames = Array.from(new Set(links.map(otherName))).sort();
    const otherChroms = otherNames.map((n) => [n, otherSizes[n]]);
    const totalOther = otherChroms.reduce((s, [, sz]) => s + sz, 0);
    const gap = Math.max(Math.round(0.01 * Math.max(totalOther, pivotSize)), 1);
    const [otherOffsets, otherSpan] = SYN.buildOffsetsLinear(otherChroms, gap);
    const xMax = Math.max(otherSpan, pivotSize);

    const barXs = [], barYs = [], barFill = [];
    const labelX = [], labelY = [], labelText = [], labelColor = [];
    for (const [name, size] of otherChroms) {
        const x0 = otherOffsets[name];
        const color = isPivotSubject ? SYN.data.targetGrey : SYN.paletteColor(SYN.data.subjectColorIndex[name], k);
        barXs.push([x0, x0 + size, x0 + size, x0]);
        barYs.push([otherRowY, otherRowY, otherRowY + barH, otherRowY + barH]);
        barFill.push(color);
        labelX.push(x0 + size / 2);
        labelY.push(otherRowY + barH / 2);
        labelText.push(name);
        labelColor.push(SYN.textColorFor(color));
        addGapTicks(otherRole, name, x0, otherRowY, xMax);
    }
    barXs.push([0, pivotSize, pivotSize, 0]);
    barYs.push([pivotRowY, pivotRowY, pivotRowY + barH, pivotRowY + barH]);
    barFill.push(pivotColor);
    labelX.push(pivotSize / 2);
    labelY.push(pivotRowY + barH / 2);
    labelText.push(pivotName);
    labelColor.push(SYN.textColorFor(pivotColor));
    addGapTicks(pivotRole, pivotName, 0, pivotRowY, xMax);

    const otherEdge = isPivotSubject ? botY + barH : topY;
    const pivotEdge = isPivotSubject ? topY : botY + barH;
    const maxScore = links.reduce((m, l) => Math.max(m, l.score), 1);
    const ribXs = [], ribYs = [], ribFill = [], ribAlpha = [], ribLabel = [];
    const sorted = links.slice().sort((a, b) => a.score - b.score);
    for (const l of sorted) {
        const oName = otherName(l);
        const ox0 = otherOffsets[oName] + otherStart(l);
        const ox1 = otherOffsets[oName] + otherEnd(l);
        const px0 = pivotStart(l), px1 = pivotEnd(l);
        ribXs.push([ox0, ox1, px1, px0]);
        ribYs.push([otherEdge, otherEdge, pivotEdge, pivotEdge]);
        const refColor = isPivotSubject ? pivotColor : SYN.paletteColor(SYN.data.subjectColorIndex[oName], k);
        ribFill.push(SYN.lighten(refColor));
        ribAlpha.push(0.25 + 0.55 * (l.score / maxScore));
        ribLabel.push(SYN.formatLinkLabel(l));
    }

    return {
        bars: {xs: barXs, ys: barYs, fill_color: barFill},
        ribbons: {xs: ribXs, ys: ribYs, fill_color: ribFill, alpha: ribAlpha, label: ribLabel},
        labels: {x: labelX, y: labelY, text: labelText, color: labelColor},
        gaps: {xs: gapXs, ys: gapYs, label: gapLabel},
        x_range: [-xMax * 0.03, xMax * 1.05],
        y_range: [botY - 0.35, topY + barH + 0.35],
        title: `${pivotName} (${pivotKind}) — ${links.length} link(s) vs ${otherNames.length} ${otherKind} chromosome(s)`,
    };
};

// One dotplot grid cell (a single target chromosome x single reference
// chromosome pair) -- always exactly two bars, no packing needed since
// there's only ever one chromosome per side. Unlike buildDetailData, this
// has to handle the zero-link case as a normal, expected outcome (not an
// edge case): the whole point of making grid cells clickable, including
// empty ones, is that "no synteny between this pair" is itself an answer
// worth showing rather than something the UI silently refuses to select.
// No orientation-based flip on ribbon x-coordinates, matching
// buildDetailData's existing convention (see its ribXs/ribYs above).
SYN.buildPairDetailData = function(targetName, subjectName, k, minScore) {
    const allLinks = SYN.data.linksByQuery[targetName] || [];
    const links = allLinks.filter((l) => l.s_chrom === subjectName && l.score >= minScore);
    const targetSize = SYN.data.querySizes[targetName];
    const subjectSize = SYN.data.subjectSizes[subjectName];
    const subjectColor = SYN.paletteColor(SYN.data.subjectColorIndex[subjectName], k);

    const barH = 0.32, topY = 1.2, botY = 0.0;
    const xMax = Math.max(targetSize, subjectSize);
    const barXs = [[0, subjectSize, subjectSize, 0], [0, targetSize, targetSize, 0]];
    const barYs = [[topY, topY, topY + barH, topY + barH], [botY, botY, botY + barH, botY + barH]];
    const barFill = [subjectColor, SYN.data.targetGrey];
    const labelX = [subjectSize / 2, targetSize / 2];
    const labelY = [topY + barH / 2, botY + barH / 2];
    const labelText = [subjectName, targetName];
    const labelColor = [SYN.textColorFor(subjectColor), SYN.textColorFor(SYN.data.targetGrey)];
    const baseX = [-xMax * 0.03, xMax * 1.05];
    const baseY = [botY - 0.35, topY + barH + 0.35];

    // both bars sit at x0=0 -- unlike buildDetailData's packed "other" row,
    // there's only ever one chromosome per side here, so no per-bar offset.
    // Same widen treatment as buildDetailData's addGapTicks -- see its own
    // comment (full bar height, no vertical inset).
    const gapXs = [], gapYs = [], gapLabel = [];
    const addGapTicks = (role, name, rowY) => {
        for (const g of SYN.gapsForChrom(role, name)) {
            const [gs, ge] = SYN.widenGapBp(g.start, g.end, xMax);
            gapXs.push([gs, ge, ge, gs]);
            gapYs.push([rowY, rowY, rowY + barH, rowY + barH]);
            gapLabel.push(SYN.formatGapLabel(name, g.start, g.end));
        }
    };
    addGapTicks('reference', subjectName, topY);
    addGapTicks('target', targetName, botY);
    const gaps = {xs: gapXs, ys: gapYs, label: gapLabel};

    if (links.length === 0) {
        return {
            bars: {xs: barXs, ys: barYs, fill_color: barFill},
            ribbons: {xs: [], ys: [], fill_color: [], alpha: [], label: []},
            labels: {x: labelX, y: labelY, text: labelText, color: labelColor},
            gaps,
            x_range: baseX, y_range: baseY,
            title: `${targetName} (${SYN.state.targetLabel}) × ${subjectName} (${SYN.state.referenceLabel}) -- no links`,
        };
    }

    const maxScore = links.reduce((m, l) => Math.max(m, l.score), 1);
    const ribXs = [], ribYs = [], ribFill = [], ribAlpha = [], ribLabel = [];
    const sorted = links.slice().sort((a, b) => a.score - b.score);
    for (const l of sorted) {
        ribXs.push([l.q_start, l.q_end, l.s_end, l.s_start]);
        ribYs.push([botY + barH, botY + barH, topY, topY]);
        ribFill.push(SYN.lighten(subjectColor));
        ribAlpha.push(0.25 + 0.55 * (l.score / maxScore));
        ribLabel.push(SYN.formatLinkLabel(l));
    }

    return {
        bars: {xs: barXs, ys: barYs, fill_color: barFill},
        ribbons: {xs: ribXs, ys: ribYs, fill_color: ribFill, alpha: ribAlpha, label: ribLabel},
        labels: {x: labelX, y: labelY, text: labelText, color: labelColor},
        gaps,
        x_range: baseX, y_range: baseY,
        title: `${targetName} (${SYN.state.targetLabel}) × ${subjectName} (${SYN.state.referenceLabel}) -- ${links.length} link(s)`,
    };
};

// mean_identity/anchor_density are null for links.tsv files predating those
// columns -- omit that line rather than printing "null". topLabel/
// bottomLabel default to the normal cross-genome case (s_chrom is always on
// the reference genome, q_chrom always on target here) -- but a homeolog
// link's two coordinates are BOTH from the same genome (see
// SYN.buildRingLayout's targetHomeologLinks/referenceHomeologLinks loops),
// so those callers pass both labels as that one genome's name instead of
// the default Reference/Target pair, which would otherwise misattribute one
// side to the wrong genome.
SYN.formatLinkLabel = function(l, topLabel, bottomLabel) {
    topLabel = topLabel || 'Reference';
    bottomLabel = bottomLabel || 'Target';
    let stats = `${l.score} genes (distinct loci)  ·  orientation=${l.orientation}`;
    if (l.mean_identity !== null && l.mean_identity !== undefined
        && l.anchor_density !== null && l.anchor_density !== undefined) {
        stats += `<br>avg identity: ${(l.mean_identity * 100).toFixed(1)}%`
            + `  ·  density: ${l.anchor_density.toFixed(1)} anchors/Mb`;
    }
    const coords = `<table style="border-collapse:collapse">`
        + `<tr><td style="padding-right:6px">${topLabel}</td>`
        + `<td style="text-align:right">${l.s_chrom}:${l.s_start.toLocaleString()}-${l.s_end.toLocaleString()}</td></tr>`
        + `<tr><td style="padding-right:6px">${bottomLabel}</td>`
        + `<td style="text-align:right">${l.q_chrom}:${l.q_start.toLocaleString()}-${l.q_end.toLocaleString()}</td></tr>`
        + `</table>`;
    return `${coords}<br>${stats}`;
};

SYN.formatGapLabel = function(chrom, start, end) {
    return `Assembly gap<br>${SYN.escapeHtml(chrom)}:${start.toLocaleString()}-${end.toLocaleString()}`
        + ` (${(end - start).toLocaleString()} bp)`;
};

// gapSource is optional (the ring/dotplot's own DoubleTap/reset handlers
// don't touch the zoom panel at all, so they never pass one) -- when given,
// writes d.gaps only if SYN.state.showGaps is on, same gating
// SYN.applyGapVisibility/SYN.applyDotplotGapVisibility use for their own
// sources, so all three panels' gap markers stay in lockstep with the one
// "Show gaps" switch.
SYN.applyDetail = function(d, barSource, ribbonSource, labelSource, detailFig, gapSource) {
    barSource.data = d.bars; barSource.change.emit();
    ribbonSource.data = d.ribbons; ribbonSource.change.emit();
    labelSource.data = d.labels; labelSource.change.emit();
    if (gapSource) {
        gapSource.data = SYN.state.showGaps ? d.gaps : {xs: [], ys: [], label: []};
        gapSource.change.emit();
    }
    detailFig.x_range.start = d.x_range[0]; detailFig.x_range.end = d.x_range[1];
    detailFig.y_range.start = d.y_range[0]; detailFig.y_range.end = d.y_range[1];
    detailFig.title.text = d.title;
    SYN.state.detailRange = {x0: d.x_range[0], x1: d.x_range[1], y0: d.y_range[0], y1: d.y_range[1]};
};

// Re-renders whatever is currently shown in the zoom panel (a pivot, a
// pair, or nothing) at a new color count / min-block-size threshold --
// shared by the color-count spinner, the palette dropdown, and the min-
// block-anchors spinner (see build_page()), so which of buildDetailData/
// buildPairDetailData applies doesn't need to be duplicated in each of
// their three CustomJS callbacks.
SYN.refreshDetail = function(k, minScore, barSource, ribbonSource, labelSource, detailFig, gapSource) {
    if (SYN.state.mode === 'pair') {
        const d = SYN.buildPairDetailData(SYN.state.pairTarget, SYN.state.pairSubject, k, minScore);
        SYN.applyDetail(d, barSource, ribbonSource, labelSource, detailFig, gapSource);
    } else if (SYN.state.pivotName) {
        const d = SYN.buildDetailData(SYN.state.pivotSide, SYN.state.pivotName, k, minScore);
        SYN.applyDetail(d, barSource, ribbonSource, labelSource, detailFig, gapSource);
    }
};

// Bokeh 3.x renders each figure inside several levels of nested Shadow DOM
// (for style encapsulation), so the actual <svg> a figure produces (with
// output_backend="svg" -- see build_page()) isn't reachable by a plain
// document.querySelector; shadow roots have to be pierced explicitly, one
// level at a time. Toolbar icons are ALSO inline <svg> elements sitting in
// the same subtree, so this collects every <svg> under the given root and
// keeps the one with the largest rendered area -- the actual plot is always
// far bigger than a ~24px icon, without needing to depend on Bokeh's
// internal class names (which aren't part of its public API and could
// change across versions).
SYN.findAllSvgsDeep = function(node, out) {
    out = out || [];
    if (!node) { return out; }
    if (node.tagName && node.tagName.toLowerCase() === 'svg') { out.push(node); return out; }
    if (node.shadowRoot) { SYN.findAllSvgsDeep(node.shadowRoot, out); }
    const children = node.children || [];
    for (let i = 0; i < children.length; i++) { SYN.findAllSvgsDeep(children[i], out); }
    return out;
};

SYN.findMainSvg = function(root) {
    const all = SYN.findAllSvgsDeep(root, []);
    let best = null, bestArea = -1;
    for (const svg of all) {
        const r = svg.getBoundingClientRect();
        const area = r.width * r.height;
        if (area > bestArea) { bestArea = area; best = svg; }
    }
    return best;
};

// Bokeh.index only holds the *root* views (here: the outer Column layout,
// plus an unrelated Notifications overlay) -- it does NOT index every nested
// view by model id, despite the name suggesting a lookup table. Reaching a
// specific Figure's view (the only way from a CustomJS callback, which only
// ever gets a plain model reference, to the actual rendered DOM/<svg>) means
// walking down from a root view through its child_views (a Map keyed by
// model, not id) until a view whose own .model.id matches turns up.
SYN.findViewByModelId = function(modelId) {
    function search(view) {
        if (!view) { return null; }
        if (view.model && view.model.id === modelId) { return view; }
        const cv = view.child_views;
        if (cv) {
            const children = typeof cv.values === 'function' ? Array.from(cv.values())
                : (Array.isArray(cv) ? cv : Object.values(cv));
            for (const child of children) {
                const found = search(child);
                if (found) { return found; }
            }
        }
        return null;
    }
    for (const root of Object.values(Bokeh.index)) {
        const found = search(root);
        if (found) { return found; }
    }
    return null;
};

// Exports a figure as a real, standalone .svg file -- no library needed at
// all, since output_backend="svg" already gives every figure a genuine
// vector <svg> DOM tree (see build_page()) rather than only ever being able
// to rasterize a canvas. This is the vector original -- open it in Inkscape,
// Illustrator, or any other vector tool to edit it or convert it to PDF. A
// rasterized alternative (PNG or JPEG, for pasting into a slide or document)
// is exported the same way via SYN.exportFigureAsRaster below, built on top
// of this same SVG -- see SYN.exportFigure for the dispatcher both save
// buttons call.
SYN.exportFigureAsSVG = function(figModel, filename) {
    const view = SYN.findViewByModelId(figModel.id);
    const svgEl = view && view.el ? SYN.findMainSvg(view.el) : null;
    if (!svgEl) {
        console.error('SYN.exportFigureAsSVG: could not find a rendered <svg> for this figure.');
        return;
    }
    // clone rather than serialize the live element directly, so the
    // namespace attribute below can be added without touching the on-page
    // original -- and so the exported file is a self-contained SVG document
    // (a bare fragment plucked out of the page has no xmlns of its own)
    const clone = svgEl.cloneNode(true);
    if (!clone.getAttribute('xmlns')) {
        clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    }
    const svgText = new XMLSerializer().serializeToString(clone);
    SYN.downloadBlob(new Blob([svgText], {type: 'image/svg+xml;charset=utf-8'}), filename);
};

// Shared by both export paths: builds a Blob URL, clicks a synthetic <a> to
// trigger the browser's download, then revokes the URL.
SYN.downloadBlob = function(blob, filename) {
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
};

// Rasterizes the same SVG SYN.exportFigureAsSVG exports, via an offscreen
// <canvas> -- no new runtime dependency, Image/Canvas are native browser
// APIs. Shared by PNG and JPEG: same pipeline, just a different mimeType
// (and, for JPEG, a quality factor -- canvas.toBlob ignores the 3rd arg for
// PNG, so passing it unconditionally is harmless). `scale` (default 2)
// renders at 2x the on-screen pixel size for a crisper-than-screen image
// rather than exposing another UI control for it. The white fill before
// drawImage matters for BOTH formats: PNG would otherwise get a transparent
// margin if Bokeh's SVG doesn't already paint an opaque background, and
// JPEG has no transparency channel at all -- without this, a transparent
// area would silently render black instead.
SYN.exportFigureAsRaster = function(figModel, filename, mimeType, scale) {
    scale = scale || 2;
    const view = SYN.findViewByModelId(figModel.id);
    const svgEl = view && view.el ? SYN.findMainSvg(view.el) : null;
    if (!svgEl) {
        console.error(`SYN.exportFigureAsRaster: could not find a rendered <svg> for this figure (${mimeType}).`);
        return;
    }
    const rect = svgEl.getBoundingClientRect();
    const clone = svgEl.cloneNode(true);
    if (!clone.getAttribute('xmlns')) {
        clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    }
    const svgText = new XMLSerializer().serializeToString(clone);
    const svgUrl = URL.createObjectURL(new Blob([svgText], {type: 'image/svg+xml;charset=utf-8'}));

    const img = new Image();
    img.onload = function() {
        const canvas = document.createElement('canvas');
        canvas.width = Math.round(rect.width * scale);
        canvas.height = Math.round(rect.height * scale);
        const ctx = canvas.getContext('2d');
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, canvas.width, canvas.height);
        ctx.scale(scale, scale);
        ctx.drawImage(img, 0, 0, rect.width, rect.height);
        URL.revokeObjectURL(svgUrl);
        canvas.toBlob(function(blob) {
            if (!blob) {
                console.error(`SYN.exportFigureAsRaster: canvas.toBlob returned null (${mimeType}).`);
                return;
            }
            SYN.downloadBlob(blob, filename);
        }, mimeType, 0.92);
    };
    img.onerror = function() {
        console.error(`SYN.exportFigureAsRaster: the rendered SVG failed to load into an <img> for rasterizing (${mimeType}).`);
        URL.revokeObjectURL(svgUrl);
    };
    img.src = svgUrl;
};

// One dispatcher both save buttons' CustomJS calls, instead of near-duplicate
// blocks each hardcoding a format.
SYN.EXPORT_MIME = {PNG: 'image/png', JPEG: 'image/jpeg'};
SYN.exportFigure = function(figModel, panel, format) {
    const filename = SYN.buildExportFilename(panel, format);
    const mimeType = SYN.EXPORT_MIME[format];
    if (mimeType) {
        SYN.exportFigureAsRaster(figModel, filename, mimeType);
    } else {
        SYN.exportFigureAsSVG(figModel, filename);
    }
};

// Builds a save button's download filename from whatever labels are
// currently active, instead of a fixed generic name, so an exported file is
// self-describing without the viewer having to rename it -- sanitized since
// the labels are free-text input and filenames can't contain arbitrary
// characters (path separators in particular). JPEG gets the conventional
// .jpg extension, not the literal (also valid, but less common) .jpeg.
SYN.buildExportFilename = function(panel, format) {
    const sanitize = (s) => (s || '').replace(/[^a-zA-Z0-9_-]+/g, '_').replace(/^_+|_+$/g, '');
    const t = sanitize(SYN.state.targetLabel) || 'target';
    const r = sanitize(SYN.state.referenceLabel) || 'reference';
    const ext = format === 'JPEG' ? 'jpg' : format.toLowerCase();
    return `${t}_vs_${r}_${panel}.${ext}`;
};

// Every title/bar-label a viewer might export reads SYN.state.targetLabel/
// referenceLabel rather than the pipeline's literal "target"/"reference"
// role tags (see SYN.buildDetailData/SYN.buildPairDetailData above) -- this
// is what actually changes those two values and refreshes the handful of
// sources/titles that don't go through a rebuild-on-every-redraw path the
// way the zoom panel does (see build_page()'s label-input callback, which
// also calls SYN.refreshDetail and SYN.applyStats after this). The 'group'
// field feeds each ruler/wedge source's hover tooltip ("Group": "@group")
// only -- nothing else reads it, so overwriting every entry is safe.
SYN.applyLabels = function(targetLabel, referenceLabel, ctx) {
    SYN.state.targetLabel = targetLabel;
    SYN.state.referenceLabel = referenceLabel;
    ctx.overview.title.text = `${referenceLabel} (reference) vs ${targetLabel} (target)`;
    ctx.dotplotFig.title.text = `${targetLabel} vs ${referenceLabel} -- click a band or a grid square to zoom`;
    for (const src of ctx.queryLikeSources) {
        src.data.group = src.data.group.map(() => targetLabel);
        src.change.emit();
    }
    for (const src of ctx.subjectLikeSources) {
        src.data.group = src.data.group.map(() => referenceLabel);
        src.change.emit();
    }
};

// The interactive HTML's replacement for the (removed) static plots' stats
// inset -- see compute_alignment_stats.py and build_page()'s --stats. Lives
// as a normal page element rather than baked into any one exportable figure,
// since export is per-panel now (ring, zoom, or dotplot), not "the whole
// image" the way a static render was.
// truncates (not rounds) to 1 decimal place -- mirrors trunc1() in this
// file's Python exactly, kept in sync deliberately
SYN.trunc1 = function(x) {
    return Math.floor(x * 10) / 10;
};

SYN.buildStatsHtml = function(targetLabel, referenceLabel) {
    const s = SYN.data.alignmentStats;
    if (!s) { return ''; }
    const pct = (n) => s.proteome_total ? SYN.trunc1(n / s.proteome_total * 100).toFixed(1) : '0.0';
    const row = (label, aligned, identity) => `
        <tr><td style="padding-right:20px">${SYN.escapeHtml(label)}</td>
            <td style="text-align:right;padding-right:20px">${aligned.toLocaleString()} (${pct(aligned)}%)</td>
            <td style="text-align:right">${SYN.trunc1(identity * 100).toFixed(1)}%</td></tr>`;
    const originLine = s.proteome_origin
        ? `<div style="color:#666;margin-top:2px">Proteome: ${SYN.escapeHtml(s.proteome_origin)}</div>` : '';
    // display:block + an explicit px width (SYN.data.statsPanelWidth, same
    // value as stats_div's own width=STATS_PANEL_WIDTH on the Python side --
    // see format_stats_html) rather than shrink-to-fit display:inline-block,
    // so this box spans the zoom panel's drawn frame (not its full nominal
    // width, which includes the toolbar strip -- see PLOT_TOOLBAR_WIDTH)
    // instead of only its own text/table content -- kept in sync with
    // format_stats_html deliberately.
    return `
        <div style="font-size:13px;color:#333;border:1px solid #ddd;border-radius:6px;
                    padding:8px 12px;display:block;box-sizing:border-box;
                    width:${SYN.data.statsPanelWidth}px;margin:4px 0 8px 0">
            <b>Alignment summary</b><br>${s.proteome_total.toLocaleString()} input proteins
            ${originLine}
            <table style="border-collapse:collapse;margin-top:4px">
                <tr style="color:#666"><th style="text-align:left;padding-right:20px"></th>
                    <th style="text-align:right;padding-right:20px">aligned</th>
                    <th style="text-align:right">avg identity</th></tr>
                ${row(targetLabel, s.query_aligned, s.query_mean_identity)}
                ${row(referenceLabel, s.subject_aligned, s.subject_mean_identity)}
            </table>
        </div>`;
};

SYN.applyStats = function(targetLabel, referenceLabel, statsDiv) {
    statsDiv.text = SYN.buildStatsHtml(targetLabel, referenceLabel);
};

// ---- client-side chaining (docs/specs/hit_table.md): decode the embedded
// hit tables, chain them (in a Web Worker when available) at load and on
// every min-identity/max-gap/hit-rank change, and hand the result to the
// same buildRingLayout/buildDotplotSegmentsForLayout/buildDetailData
// functions above -- no separate render path for "the first chain" vs
// "a re-chain after a control change". ----

// the chainer is always run at this minBlock (docs/specs/hit_table.md
// section 4.6: extraction doesn't depend on it, so this is the least
// strict floor that still excludes single/double-anchor noise) -- the min
// block size control is a pure client-side filter on top of the reply (see
// SYN.data.ribbons' own comment above), never a reason to re-chain.
SYN.CHAIN_REQUEST_MIN_BLOCK = 3;

SYN.HIT_RANK_MAP = {'best only': 1, '≤ 2': 2, '≤ 3': 3, all: SYNCHAIN.DEFAULTS.maxHitRank};

SYN.now = function() {
    return (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
};

// Reads the three chain-affecting controls (min identity/max gap/hit rank)
// off SYN.ui -- the min block size control is deliberately NOT here (see
// SYN.CHAIN_REQUEST_MIN_BLOCK above): it never changes what gets requested,
// only what SYN.applyChainResult's callers filter the reply down to.
SYN.currentChainParams = function() {
    return {
        minPositive: SYN.ui.minIdentitySpinner.value / 100,
        maxGap: SYN.ui.maxGapSpinner.value,
        maxHitRank: SYN.HIT_RANK_MAP[SYN.ui.hitRankSelect.value],
        minBlock: SYN.CHAIN_REQUEST_MIN_BLOCK,
    };
};

// Runs inside the Worker, prepended with the inlined chain.js source (see
// SYN.startChainer) and NOTHING else from this file -- a Worker gets its own
// fresh global scope, with none of SHARED_JS's own SYN.* functions in it, so
// this calls SYNCHAIN directly rather than through SYN.buildChainers/
// SYN.runChainers below (the main-thread fallback's own equivalent, which
// CAN use them since it runs in this same scope). The ${...} substitution
// happens once, in SHARED_JS itself, when this string is built -- by the
// time it reaches a Worker it's a plain literal, not a template needing
// anything of SYN's still in scope. Kept as one string (rather than a named
// function serialized via toString()) so it reads the same whichever path
// runs it, and so a syntax error in it fails at page-build time (Python
// would embed broken JS) rather than silently inside a Worker.
SYN.WORKER_SHIM_SRC = `
let chainers = null;
self.onmessage = function(e) {
    const msg = e.data;
    if (msg.type === 'init') {
        chainers = {
            cross: SYNCHAIN.createChainer(msg.target, msg.reference, {}),
            targetSelf: SYNCHAIN.createChainer(msg.target, msg.target, {selfMode: true}),
            referenceSelf: SYNCHAIN.createChainer(msg.reference, msg.reference, {selfMode: true}),
        };
        return;
    }
    if (msg.type === 'run') {
        const p = Object.assign({}, msg.params, {minBlock: ${SYN.CHAIN_REQUEST_MIN_BLOCK}});
        const t0 = performance.now();
        const cross = chainers.cross.run(p);
        const targetSelf = chainers.targetSelf.run(p);
        const referenceSelf = chainers.referenceSelf.run(p);
        self.postMessage({
            seq: msg.seq,
            crossLinks: SYNCHAIN.blocksToLinks(cross),
            targetSelfLinks: SYNCHAIN.blocksToLinks(targetSelf),
            referenceSelfLinks: SYNCHAIN.blocksToLinks(referenceSelf),
            ms: performance.now() - t0,
        });
    }
};
`;

// Three chainers sharing the decoded tables: cross (target vs reference)
// plus each genome's own self-comparison (homeologs) -- createChainer's own
// per-stage caching (docs/specs/hit_table.md section 3's own note on this)
// means a min-identity/max-gap/hit-rank change that doesn't move the anchor
// set at all (rare, but e.g. a repeated hit-rank click) costs next to
// nothing on a later .run() call. Used only by SYN.startChainer's
// main-thread fallback -- the Worker path above inlines the equivalent
// logic itself, since it has no SYN.* of its own to call into.
SYN.buildChainers = function(target, reference) {
    return {
        cross: SYNCHAIN.createChainer(target, reference, {}),
        targetSelf: SYNCHAIN.createChainer(target, target, {selfMode: true}),
        referenceSelf: SYNCHAIN.createChainer(reference, reference, {selfMode: true}),
    };
};

SYN.runChainers = function(chainers, params, seq) {
    const p = Object.assign({}, params, {minBlock: SYN.CHAIN_REQUEST_MIN_BLOCK});
    const t0 = SYN.now();
    const cross = chainers.cross.run(p);
    const targetSelf = chainers.targetSelf.run(p);
    const referenceSelf = chainers.referenceSelf.run(p);
    return {
        seq,
        crossLinks: SYNCHAIN.blocksToLinks(cross),
        targetSelfLinks: SYNCHAIN.blocksToLinks(targetSelf),
        referenceSelfLinks: SYNCHAIN.blocksToLinks(referenceSelf),
        ms: SYN.now() - t0,
    };
};

// Decodes the embedded payload and sets up however this browser will chain
// it -- a Web Worker (built from the SAME inlined chain.js source text the
// main thread already has, read back out of the inert <script
// type="text/plain" id="synchain-src"> copy build_page() emits, plus
// SYN.WORKER_SHIM_SRC) so a re-chain never blocks the UI thread, or, if
// constructing the Worker throws (e.g. a `data:`/`file:` page in a browser
// that restricts Worker creation there -- see the plan's own note), a
// synchronous main-thread fallback with exactly one console.warn, never a
// hard failure.
SYN.startChainer = async function(payload) {
    const chroms = {
        target: {names: SYN.data.queryNamesNatural, sizes: SYN.data.queryNamesNatural.map((n) => SYN.data.querySizes[n])},
        reference: {names: SYN.data.subjectNamesNatural, sizes: SYN.data.subjectNamesNatural.map((n) => SYN.data.subjectSizes[n])},
    };
    const tables = await SYNCHAIN.decodePayload(payload, chroms);
    SYN.chain = {
        tables, seq: -1, appliedSeq: -1, busy: false,
        pendingSeq: null, pendingParams: null, worker: null, useWorker: true, mainThreadChainers: null,
    };
    try {
        const chainSrc = document.getElementById('synchain-src').textContent;
        const blob = new Blob([chainSrc, SYN.WORKER_SHIM_SRC], {type: 'text/javascript'});
        SYN.chain.worker = new Worker(URL.createObjectURL(blob));
        SYN.chain.worker.onmessage = function(e) { SYN.onChainMessage(e.data); };
        SYN.chain.worker.postMessage({type: 'init', target: tables.target, reference: tables.reference});
    } catch (e) {
        console.warn('SYN.startChainer: Worker unavailable, chaining on the main thread instead.', e);
        SYN.chain.useWorker = false;
        SYN.chain.mainThreadChainers = SYN.buildChainers(tables.target, tables.reference);
    }
};

// While a request is pending and nothing is zoomed, the zoom panel's own
// title doubles as a busy indicator -- restored to its normal (possibly
// still-empty) state by SYN.applyChainResult once the reply lands.
SYN.showComputingStatus = function() {
    if (SYN.state.mode === null && !SYN.state.pivotName) {
        SYN.ui.detailFig.title.text = 'Computing synteny…';
    }
};

// At most one chain request in flight at a time -- while busy, only the
// latest params are remembered (not queued), and sent the moment the
// current reply arrives (see SYN.onChainMessage). seq is used to ignore any
// reply older than the latest one actually applied, which can otherwise
// happen if the main-thread fallback and a slow Worker reply somehow race
// (they can't in practice -- one browser uses exactly one path -- but the
// check is free and matches the plan's own wording).
SYN.requestChain = function(params) {
    SYN.chain.seq++;
    const seq = SYN.chain.seq;
    if (SYN.chain.busy) {
        SYN.chain.pendingSeq = seq;
        SYN.chain.pendingParams = params;
        return;
    }
    SYN.chain.busy = true;
    SYN.showComputingStatus();
    if (SYN.chain.useWorker) {
        SYN.chain.worker.postMessage({type: 'run', seq, params});
    } else {
        SYN.onChainMessage(SYN.runChainers(SYN.chain.mainThreadChainers, params, seq));
    }
};

SYN.onChainMessage = function(msg) {
    SYN.chain.busy = false;
    if (msg.seq > SYN.chain.appliedSeq) {
        SYN.chain.appliedSeq = msg.seq;
        SYN.applyChainResult(msg);
    }
    if (SYN.chain.pendingSeq !== null) {
        const seq = SYN.chain.pendingSeq, params = SYN.chain.pendingParams;
        SYN.chain.pendingSeq = null;
        SYN.chain.pendingParams = null;
        SYN.chain.busy = true;
        SYN.showComputingStatus();
        if (SYN.chain.useWorker) {
            SYN.chain.worker.postMessage({type: 'run', seq, params});
        } else {
            SYN.onChainMessage(SYN.runChainers(SYN.chain.mainThreadChainers, params, seq));
        }
    }
};

SYN.groupLinksBy = function(links, key) {
    const g = {};
    for (const l of links) {
        (g[l[key]] || (g[l[key]] = [])).push(l);
    }
    return g;
};

// Rebuilds SYN.data's link tables from a fresh chain result and re-renders
// every panel through the exact same apply* functions every other control
// change already uses (SYN.applyRingLayout, SYN.applyDotplotOrder/
// SYN.applyDotplotSegmentsForCurrentLayout, SYN.refreshDetail) -- this is
// the one place a chain result ever reaches SYN.data, whether it came from
// the very first load or the Nth min-identity change.
SYN.applyChainResult = function(msg) {
    const ui = SYN.ui;
    // kept flat (sorted, chain.js's own compareBlocks order -- see
    // docs/specs/hit_table.md section 5) so SYN.exportBlocksTsv can filter
    // and re-serialize it without needing to flatten linksByQuery back out
    // in some arbitrary (and here, non-deterministic: Object.values order)
    // order of its own.
    SYN.data.crossLinksFlat = msg.crossLinks;
    SYN.data.linksByQuery = SYN.groupLinksBy(msg.crossLinks, 'q_chrom');
    SYN.data.linksByReference = SYN.groupLinksBy(msg.crossLinks, 's_chrom');
    SYN.data.targetHomeologLinks = msg.targetSelfLinks;
    SYN.data.referenceHomeologLinks = msg.referenceSelfLinks;
    let maxScore = 1;
    for (const l of msg.crossLinks) { if (l.score > maxScore) { maxScore = l.score; } }
    for (const l of msg.targetSelfLinks) { if (l.score > maxScore) { maxScore = l.score; } }
    for (const l of msg.referenceSelfLinks) { if (l.score > maxScore) { maxScore = l.score; } }
    SYN.data.maxScore = maxScore;

    const ringOrder = SYN.ringOrderFor(SYN.state.orderBySize);
    SYN.data.ribbons = SYN.buildRibbonRecords(SYN.state.ringQueryOffsets, SYN.state.ringSubjectOffsets);
    SYN.applyOverviewRibbons(ui.minBlockSpinner.value, ui.colorSpinner.value, ui.selfLinksToggle.active,
                              ui.hideSyntenyToggle.active, ui.ribbonSource);
    if (SYN.state.orderBySimilarity) {
        const dpOrder = SYN.dpOrderFor(true);
        SYN.applyDotplotOrder(dpOrder.queryOrder, dpOrder.subjectOrder, ui.colorSpinner.value, ui.minBlockSpinner.value, {
            query: ui.dpQuerySource, subject: ui.dpSubjectSource, grid: ui.dpGridSource,
            queryLabel: ui.dpQueryLabelSource, subjectLabel: ui.dpSubjectLabelSource,
            segment: ui.dpSegmentSource, gap: ui.dpGapSource, fig: ui.dotplotFig,
        });
    } else {
        SYN.applyDotplotSegmentsForCurrentLayout(ui.minBlockSpinner.value, ui.colorSpinner.value, ui.dpSegmentSource);
    }
    SYN.refreshDetail(ui.colorSpinner.value, ui.minBlockSpinner.value, ui.detailBarSource, ui.detailRibbonSource,
                       ui.detailLabelSource, ui.detailFig, ui.detailGapSource);
    // SYN.refreshDetail is a no-op with nothing zoomed -- restore the panel
    // title SYN.showComputingStatus overwrote for the duration of this
    // request, same condition that function itself gates on.
    if (SYN.state.mode === null && !SYN.state.pivotName) {
        ui.detailFig.title.text = SYN.emptyDetail().title;
    }
    ui.minBlockSpinner.high = maxScore;

    SYN.chain.lastMs = msg.ms;
    SYN.updateChainStatus();
};

// "N block(s) · X ms": blocks currently passing the min-block filter, and how
// long the last chain request took -- called after every chain result AND
// every min-block change (a pure filter that never re-chains).
SYN.updateChainStatus = function() {
    const ui = SYN.ui;
    if (!SYN.data.crossLinksFlat || !SYN.chain) { return; }
    const visible = SYN.data.crossLinksFlat.filter((l) => l.score >= ui.minBlockSpinner.value).length;
    ui.chainStatusDiv.text = `${visible} block(s) · ${(SYN.chain.lastMs || 0).toFixed(0)} ms`;
};

// "the cross blocks currently on screen (min block size filter applied)"
// (the plan's own wording) -- SYN.data.crossLinksFlat is already sorted in
// exactly the order docs/specs/hit_table.md section 5 specifies (it's
// chain.js's own blocksToLinks output, untouched -- see SYN.applyChainResult),
// and filtering a sorted list by a monotonic predicate preserves that order,
// so this is byte-identical to a fresh `node bin/chain_blocks.mjs --min_block
// <current value>` run with the same min-identity/max-gap/hit-rank.
SYN.exportBlocksTsv = function() {
    const minScore = SYN.ui.minBlockSpinner.value;
    const filtered = (SYN.data.crossLinksFlat || []).filter((l) => l.score >= minScore);
    const tsv = SYNCHAIN.linksToTsv(filtered);
    SYN.downloadBlob(new Blob([tsv], {type: 'text/tab-separated-values'}), SYN.buildBlocksExportFilename());
};

SYN.buildBlocksExportFilename = function() {
    const sanitize = (s) => (s || '').replace(/[^a-zA-Z0-9_-]+/g, '_').replace(/^_+|_+$/g, '');
    const t = sanitize(SYN.state.targetLabel) || 'target';
    const r = sanitize(SYN.state.referenceLabel) || 'reference';
    const id = Math.round(SYN.ui.minIdentitySpinner.value);
    return `${t}_vs_${r}_blocks_id${id}_gap${SYN.ui.maxGapSpinner.value}_min${SYN.ui.minBlockSpinner.value}.tsv`;
};

// The page's one entry point, run once from build_page()'s
// doc.js_on_event(DocumentReady, ...) -- see that call for which
// sources/figure/widget end up on `s`. Stashes `s` as SYN.ui (the registry
// every later chain-related function above reads/writes through), does an
// initial ring/dotplot render with zero links (through the exact same
// functions every later control change reuses -- there is no separate
// "initial render" code path, the way Python's direct geometry construction
// used to be -- see module docstring), then kicks off the actual chaining:
// decode the embedded hit tables, seed the min-identity/min-block controls
// from --min_identity/--min_block if given or SYNCHAIN.autoParams otherwise
// (again, before a viewer could possibly have touched either), and request
// the first real chain.
SYN.init = function(s) {
    SYN.ui = s;
    SYN.state.targetLabel = SYN.data.targetLabelDefault;
    SYN.state.referenceLabel = SYN.data.referenceLabelDefault;

    const ringOrder = SYN.ringOrderFor(true);
    SYN.applyRingLayout(ringOrder.queryOrder, ringOrder.subjectOrder, {
        querySource: s.querySource, subjectSource: s.subjectSource, labelSource: s.labelSource,
        ribbonSource: s.ribbonSource, gapSource: s.gapSource,
        minScore: s.minBlockSpinner.value, k: s.colorSpinner.value,
        showSelfLinks: false, hideSynteny: false,
    });

    const dpOrder = SYN.dpOrderFor(false);
    SYN.applyDotplotOrder(dpOrder.queryOrder, dpOrder.subjectOrder, s.colorSpinner.value, s.minBlockSpinner.value, {
        query: s.dpQuerySource, subject: s.dpSubjectSource, grid: s.dpGridSource,
        queryLabel: s.dpQueryLabelSource, subjectLabel: s.dpSubjectLabelSource,
        segment: s.dpSegmentSource, gap: s.dpGapSource, fig: s.dotplotFig,
    });

    SYN.startChainer(SYN.hitsPayload).then(function() {
        const auto = SYNCHAIN.autoParams(SYN.chain.tables.target, SYN.chain.tables.reference);
        const init = SYN.data.initialParams;
        const minPositive = init.minPositive != null ? init.minPositive : auto.minPositive;
        const minBlock = init.minBlock != null ? init.minBlock : auto.minBlock;
        // setting these fires min_identity_spinner/min_block_spinner's own
        // change callbacks (see build_page()) -- harmless here: SYN.data
        // still holds the empty placeholders the ring/dotplot render above
        // already rendered from, and the real SYN.requestChain call right
        // after this makes both redundant re-renders moot within a frame.
        s.minIdentitySpinner.value = Math.round(minPositive * 100);
        s.minBlockSpinner.value = minBlock;
        SYN.requestChain(SYN.currentChainParams());
    }).catch(function(e) {
        console.error('SYN.startChainer failed -- the page has no synteny to show:', e);
    });
};

"""


def build_page(ds, query_name, subject_name, query_subtitle=None, subject_subtitle=None,
               alignment_stats=None, hits_payload=None, chain_js_src='',
               min_identity=None, max_gap=25, min_block=None):
    query_names = [n for n, _ in ds.query_chroms]
    subject_names = [n for n, _ in ds.subject_chroms]
    subject_index = {name: i for i, (name, _) in enumerate(ds.subject_chroms)}
    query_index = {name: i for i, (name, _) in enumerate(ds.query_chroms)}

    q_src, s_src, r_src, label_src = build_overview_sources()

    # Assembly-gap wedges (both genomes' own gaps, drawn on their own half of
    # the ring): built once, client-side, by SYN.buildRingLayout at document-
    # ready (see SYN.init below) into SYN.data.gapRecords, replaced wholesale
    # on every reorder after that. gap_src itself starts empty --
    # SYN.applyGapVisibility (also run from SYN.init) is what actually copies
    # gapRecords into it, gated on show_gaps_toggle's own default (False, see
    # below), so the two can never disagree the way a Python-rendered initial
    # state could.
    gap_src = ColumnDataSource(dict(xs=[], ys=[], label=[]))

    lim = OUTER_R + 0.22
    # no 'tap' here -- that shorthand adds its own generic, unrestricted
    # TapTool, which would sit in the toolbar (visible) alongside the
    # renderer-scoped, hidden one added below via add_tools(), showing two
    # redundant Tap buttons for the same gesture. The renderer-scoped one
    # is the one actually driving click-to-zoom (see .selected.js_on_change
    # wiring below) and is the only one that needs to exist at all.
    overview = figure(width=620, height=620, match_aspect=True,
                       x_range=Range1d(-lim, lim), y_range=Range1d(-lim, lim),
                       title=f"{subject_subtitle or subject_name} (reference) vs "
                             f"{query_subtitle or query_name} (target)",
                       tools="pan,wheel_zoom,reset",
                       output_backend="svg")
    overview.axis.visible = False
    overview.grid.visible = False
    overview.toolbar.logo = None

    ribbon_renderer = overview.patches('xs', 'ys', source=r_src, fill_color='fill_color',
                                        line_color=None, fill_alpha='alpha')
    query_renderer = overview.patches('xs', 'ys', source=q_src, fill_color='fill_color',
                                       line_color='black', line_width=0.5)
    subject_renderer = overview.patches('xs', 'ys', source=s_src, fill_color='fill_color',
                                         line_color='black', line_width=0.5)
    overview.text('x', 'y', source=label_src, text_align='center', text_baseline='middle',
                  text_font_size='10px')
    # drawn last (on top of the wedges) so a gap stripe is actually visible
    # against the chromosome band it marks. Not part of overview_tap below --
    # a gap marker isn't a clickable region (there's nothing to zoom into
    # that the chromosome wedge underneath doesn't already offer), just a
    # hover target -- and Bokeh hit-tests each TapTool renderer independently
    # of paint order, so leaving it out doesn't block clicks to the wedge
    # underneath either.
    gap_renderer = overview.patches('xs', 'ys', source=gap_src, fill_color='#000000',
                                     fill_alpha=0.55, line_color=None)

    # follow_mouse -- see detail_fig's identical HoverTool below for why: a
    # ribbon can span most of the ring, so the default snap-to-center
    # tooltip can land outside the current (zoomed/panned) view
    overview.add_tools(HoverTool(renderers=[ribbon_renderer], tooltips="@label{safe}",
                                  point_policy='follow_mouse'))
    overview.add_tools(HoverTool(renderers=[query_renderer, subject_renderer],
                                  tooltips=[("Chromosome", "@name"), ("Size", "@size_label"),
                                            ("Group", "@group")]))
    overview.add_tools(HoverTool(renderers=[gap_renderer], tooltips="@label{safe}",
                                  point_policy='follow_mouse'))
    # both sides are clickable: either wedge triggers the same zoom behavior.
    # ribbon_renderer is on the same TapTool (not a second one) -- clicking a
    # ribbon highlights it instead of zooming (see SYN.applyRibbonHighlight
    # and its wiring below); Bokeh hit-tests per-renderer within one TapTool
    # just fine, and a second TapTool would mean a second toolbar button (see
    # the "no 'tap' here" note above -- same reasoning against duplicates).
    # visible=False hides its toolbar button -- clicking a wedge doesn't
    # need one, it's not a mode you toggle on/off -- but a tool added via
    # add_tools() (unlike one named in the tools= string above) is never
    # auto-activated, hidden or not, so active_tap has to be set explicitly
    # or the hidden tool would sit there never actually receiving clicks
    overview_tap = TapTool(renderers=[subject_renderer, query_renderer, ribbon_renderer], visible=False)
    overview.add_tools(overview_tap)
    overview.toolbar.active_tap = overview_tap
    # double-click to reset pan/zoom -- a shortcut for the toolbar's own
    # Reset button, since this ring never reprograms its own range after
    # creation, the reset target is just its own original (fixed) extent.
    # Also clears an active ribbon highlight (see SYN.applyRibbonHighlight) --
    # a double-click reads as "reset this panel", which should include any
    # transient ribbon selection, not just pan/zoom.
    overview.js_on_event(DoubleTap, CustomJS(args=dict(fig=overview, ribbon_source=r_src), code=f"""
        fig.x_range.start = {-lim}; fig.x_range.end = {lim};
        fig.y_range.start = {-lim}; fig.y_range.end = {lim};
        ribbon_source.selected.indices = [];
    """))

    # Whole-genome dotplot: target along x, reference along y, each block
    # drawn as the line segment its own (query, subject) span traces out --
    # the standard synteny-dotplot convention (a diagonal streak per
    # collinear block, anti-diagonal for an inversion). The chromosome-ruler
    # strips along each axis are clickable, wired into the exact same
    # zoom-pivot mechanism the ring wedges use below (same 'name' field, same
    # tap callback) -- from SYN.buildDetailData's perspective, a click is a
    # click regardless of which panel it came from.
    (dp_q_src, dp_s_src, dp_seg_src, dp_grid_src, dp_q_label_src, dp_s_label_src) = build_dotplot_sources()

    # Dotplot gap lines (nice-to-have companion to the ring's gap wedges,
    # same "Show gaps" switch) -- entirely client-side, same as the segments
    # themselves: SYN.buildDotplotGapLines rebuilds it from
    # SYN.data.targetGapsByChrom/referenceGapsByChrom + whichever dotplot
    # order is currently active whenever it's needed. Starts empty, matching
    # the switch's default-off state.
    dp_gap_src = ColumnDataSource(dict(xs=[], ys=[], label=[]))

    # Placeholder extent -- SYN.init's first call to SYN.applyDotplotOrder
    # (see build_page()'s doc.js_on_event(DocumentReady, ...) below) sets the
    # real x_range/y_range before the page ever paints, the same way every
    # later reorder/min-length-filter change does (see that function's own
    # comment on the *3/*1.02 multipliers).
    # no 'tap' here -- see overview's identical figure(...) above for why
    dotplot_fig = figure(width=620, height=620,
                          x_range=Range1d(-1, 1), y_range=Range1d(-1, 1),
                          title="Whole-genome dotplot -- click a band or a grid square to zoom",
                          tools="pan,wheel_zoom,reset",
                          output_backend="svg")
    dotplot_fig.axis.visible = False
    dotplot_fig.grid.visible = False
    dotplot_fig.toolbar.logo = None

    dotplot_fig.multi_line('xs', 'ys', source=dp_grid_src, line_color='black', line_width=1)
    dp_segment_renderer = dotplot_fig.multi_line('xs', 'ys', source=dp_seg_src, line_color='line_color',
                                                  line_alpha='alpha', line_width=2)
    # thin, light-grey dotted lines, drawn above the segments/gridlines but
    # below the ruler strips -- a gap line crossing the ruler itself would
    # look like noise against the strip's own solid fill, not add
    # information there. Deliberately subtle (thinner and lighter than the
    # ring/zoom panel's own black gap markers): the dotplot is already the
    # busiest panel, and a gap line here is a secondary "here's roughly
    # where" cue, not the primary way to inspect one -- see GAP_LINE_COLOR's
    # own comment.
    dp_gap_renderer = dotplot_fig.multi_line('xs', 'ys', source=dp_gap_src, line_color=GAP_LINE_COLOR,
                                              line_dash='dotted', line_width=0.75)
    dp_query_renderer = dotplot_fig.patches('xs', 'ys', source=dp_q_src, fill_color='fill_color',
                                             line_color='black', line_width=0.5)
    dp_subject_renderer = dotplot_fig.patches('xs', 'ys', source=dp_s_src, fill_color='fill_color',
                                               line_color='black', line_width=0.5)
    dotplot_fig.text('x', 'y', source=dp_q_label_src, text_align='center', text_baseline='top',
                      text_font_size='8px')
    dotplot_fig.text('x', 'y', source=dp_s_label_src, text_align='right', text_baseline='middle',
                      text_font_size='8px')

    # follow_mouse -- see overview's identical HoverTool above for why
    dotplot_fig.add_tools(HoverTool(renderers=[dp_segment_renderer], tooltips="@label{safe}",
                                     point_policy='follow_mouse'))
    dotplot_fig.add_tools(HoverTool(renderers=[dp_query_renderer, dp_subject_renderer],
                                     tooltips=[("Chromosome", "@name"), ("Size", "@size_label"),
                                               ("Group", "@group")]))
    # line_policy='interp' (default is 'nearest', which -- since each gap
    # line here is just 2 vertices, its own start and end -- would otherwise
    # snap the tooltip to whichever of those two endpoints is closer, so it
    # jumps between the top/bottom (or left/right) of the plot instead of
    # tracking the mouse along the rest of the line's length; interp gives a
    # continuously-interpolated position along the segment for
    # point_policy='follow_mouse' to actually follow)
    dotplot_fig.add_tools(HoverTool(renderers=[dp_gap_renderer], tooltips="@label{safe}",
                                     point_policy='follow_mouse', line_policy='interp'))
    # one TapTool covering both ruler renderers (rather than one per
    # renderer, as an earlier version had) -- a single tool can span
    # renderers backed by different data sources just fine, each click still
    # only touches whichever renderer/glyph was actually hit, and it means
    # only one tool's active_tap needs setting below. A grid-square tap (the
    # pairwise zoom) is handled separately, below, by a plain Tap event
    # handler and SYN.dpCellAt -- there is no invisible per-cell renderer
    # here any more for a TapTool to hit-test against (see module docstring).
    # visible=False -- see overview's identical TapTool above for why.
    dotplot_tap = TapTool(renderers=[dp_query_renderer, dp_subject_renderer], visible=False)
    dotplot_fig.add_tools(dotplot_tap)
    dotplot_fig.toolbar.active_tap = dotplot_tap
    # double-click to reset -- see overview's identical handler above.
    # Reads SYN.state's CURRENT rx/ry/totalX/totalY (set for real by
    # SYN.init, see this function's own doc.js_on_event(DocumentReady, ...)
    # below) rather than a value fixed once at render time: a reorder alone
    # never changes them, but the min-length filter (see SYN.filterBySize)
    # can shrink them, and resetting to the ORIGINAL (pre-filter) span would
    # leave dead margin beyond wherever the plot now actually ends.
    dotplot_fig.js_on_event(DoubleTap, CustomJS(args=dict(fig=dotplot_fig), code="""
        fig.x_range.start = -SYN.state.dpRx * 3; fig.x_range.end = SYN.state.dpTotalX * 1.02;
        fig.y_range.start = -SYN.state.dpRy * 3; fig.y_range.end = SYN.state.dpTotalY * 1.02;
    """))

    detail_bar_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[]))
    detail_rib_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[], alpha=[], label=[]))
    detail_label_src = ColumnDataSource(dict(x=[], y=[], text=[], color=[]))
    detail_gap_src = ColumnDataSource(dict(xs=[], ys=[], label=[]))
    detail_fig = figure(width=DETAIL_FIG_WIDTH, height=302, x_axis_label='position (Mb)',
                         title="Click any chromosome wedge on the left to zoom in.",
                         tools="pan,wheel_zoom,reset",
                         output_backend="svg")
    detail_fig.yaxis.visible = False
    detail_fig.toolbar.logo = None
    # the underlying bar/ribbon coordinates stay in bp (same units build_page's
    # other geometry uses) -- only the tick labels are rescaled for display,
    # so nothing about the actual layout math needs to change
    detail_fig.xaxis.formatter = CustomJSTickFormatter(code="""
        const mb = tick / 1e6;
        return (Math.round(mb * 10) / 10).toString();
    """)
    detail_ribbon_renderer = detail_fig.patches('xs', 'ys', source=detail_rib_src,
                                                 fill_color='fill_color', line_color=None,
                                                 fill_alpha='alpha')
    detail_bar_renderer = detail_fig.patches('xs', 'ys', source=detail_bar_src,
                                              fill_color='fill_color', line_color='black', line_width=0.5)
    # drawn after the bars so a gap tick is visible on top of the bar it
    # marks -- see overview's identical gap_renderer above for why it's not
    # part of any TapTool (nothing to zoom into beyond the bar underneath)
    detail_gap_renderer = detail_fig.patches('xs', 'ys', source=detail_gap_src, fill_color='#000000',
                                              fill_alpha=0.55, line_color=None)
    detail_fig.text('x', 'y', source=detail_label_src, text_align='center', text_baseline='middle',
                     text_font_size='9px', text_color='color')
    # point_policy='follow_mouse' (default is 'snap_to_data', anchored to
    # the ribbon's geometric center) -- a ribbon can span the whole panel,
    # so anchoring to its center means the tooltip can render outside the
    # current view once someone's zoomed into one end of it
    detail_fig.add_tools(HoverTool(renderers=[detail_ribbon_renderer], tooltips="@label{safe}",
                                    point_policy='follow_mouse'))
    detail_fig.add_tools(HoverTool(renderers=[detail_gap_renderer], tooltips="@label{safe}",
                                    point_policy='follow_mouse'))
    # double-click to reset -- unlike the ring/dotplot above, this panel's
    # range is rewritten on every click (see SYN.applyDetail), so the reset
    # target has to be whatever SYN.applyDetail most recently stashed in
    # SYN.state.detailRange rather than a value fixed at page load; a no-op
    # before the first click, when there's nothing to reset to yet
    detail_fig.js_on_event(DoubleTap, CustomJS(args=dict(fig=detail_fig), code="""
        const r = SYN.state.detailRange;
        if (!r) { return; }
        fig.x_range.start = r.x0; fig.x_range.end = r.x1;
        fig.y_range.start = r.y0; fig.y_range.end = r.y1;
    """))

    # default to the actual input file name/accession (query_subtitle/
    # subject_subtitle -- see --query_subtitle/--subject_subtitle) rather
    # than the literal "target"/"reference" role tags, so a viewer sees
    # which physical genome is which immediately, right in the box itself --
    # there used to be a separate "target: <file>" Div below the header for
    # this; removed now that it would just repeat what this box already
    # says. Falls back to the role tag itself when no subtitle was given
    # (--query_subtitle/--subject_subtitle are both optional). A viewer can
    # still retype either box to anything else before exporting a panel,
    # see SYN.applyLabels.
    target_label_input = TextInput(title="Target label", value=query_subtitle or query_name,
                                    width=TOP_CONTROL_WIDTH)
    # "reference", not subject_name -- subject_name is "comparison", the
    # internal pipeline role tag (channel/filename convention, see main.nf),
    # not the user-facing term this box's own title already uses
    reference_label_input = TextInput(title="Reference label", value=subject_subtitle or "reference",
                                       width=TOP_CONTROL_WIDTH)

    palette_select = Select(title="Color palette", value=DEFAULT_PALETTE_NAME,
                             options=list(PALETTES.keys()), width=TOP_CONTROL_WIDTH)
    color_spinner = Spinner(title="Colors", low=1, high=MAX_COLORS,
                             step=1, value=DEFAULT_COLORS, width=TOP_CONTROL_WIDTH)
    # Chain parameters (docs/specs/hit_table.md section 3) -- min identity,
    # max gap, and hit rank each trigger a fresh client-side re-chain
    # (SYN.requestChain); min block size is a pure post-filter (see
    # SYN.data.ribbons' own comment above) and never re-chains. min_identity/
    # min_block's real initial values depend on the actual hit tables
    # (SYNCHAIN.autoParams), which Python never sees -- SYN.init sets both
    # from --min_identity/--min_block if given, else from autoParams, the
    # moment the embedded hit tables are decoded (before a viewer could
    # touch either control, same "no separate initial render" reasoning as
    # every other control on this page). The placeholder values below are
    # only ever visible for the instant between DOM construction and that
    # first SYN.init pass.
    min_identity_spinner = Spinner(title="Min identity (%)", low=30, high=100, step=1,
                                    value=round((min_identity or 0.5) * 100), width=TOP_CONTROL_WIDTH)
    max_gap_spinner = Spinner(title="Max gap (genes)", low=1, high=100, step=1,
                               value=max_gap, width=TOP_CONTROL_WIDTH)
    HIT_RANK_OPTIONS = ["best only", "≤ 2", "≤ 3", "all"]
    hit_rank_select = Select(title="Hit rank", value="all", options=HIT_RANK_OPTIONS, width=TOP_CONTROL_WIDTH)
    # low=3: docs/specs/hit_table.md's MIN_CHAIN_LENGTH is 2, but a 2-anchor
    # chain is barely evidence of anything -- 3 is this control's own floor,
    # independent of the (always looser) minBlock=3 every chain request
    # itself is run at (see SYN.data.ribbons' own comment above). high is a
    # generous placeholder -- SYN.applyChainResult updates it to the actual
    # largest block on every result, same "don't cap at an arbitrary round
    # number" reasoning min_seq_size_spinner's high uses below.
    min_block_spinner = Spinner(title="Min block size", low=3, high=1000,
                                 step=1, value=min_block or 5, width=TOP_CONTROL_WIDTH)
    # One-line status ("N block(s) · X ms"), refreshed by every
    # SYN.applyChainResult -- lets a viewer tell a slow re-chain (a large
    # genome, a loose max-gap) apart from "nothing matched".
    chain_status_div = Div(text="", width=TOP_CONTROL_WIDTH, margin=(18, 0, 0, 0))
    # Post-hoc chromosome-length filter for the ring + dotplot, independent
    # of --min_seq_size: that pipeline flag already dropped anything shorter
    # than its own threshold before this script ever saw the data (see
    # RENAME_SEQUENCES's chrom_sizes output, main.nf), so this control can
    # only ever filter what actually made it into the page, never recover
    # what didn't. Starts at 0 (Mb, like the ruler size labels elsewhere on
    # this page) -- 0 means "show every chromosome embedded in this page",
    # not "no filtering ever happened upstream". high is this dataset's own
    # largest chromosome -- `or 1` guards the degenerate empty-dataset case,
    # where Spinner would otherwise get low == high == 0.
    max_seq_size = max(list(ds.query_sizes.values()) + list(ds.subject_sizes.values()), default=0)
    min_seq_size_spinner = Spinner(title="Min sequence length (Mb)", low=0,
                                    high=round(max_seq_size / 1e6, 2) or 1,
                                    step=0.1, value=0, width=TOP_CONTROL_WIDTH)
    reset_btn = Button(label="✕ clear zoom", button_type="default", width=140,
                        height=TOOLBAR_CONTROL_HEIGHT)
    # label is just "save" (not "save ring"/"save zoom"/"save dotplot") and
    # the button narrow to match -- which panel it saves is already obvious
    # from its position directly under that panel, so the longer label was
    # only ever spending width, not clarity
    SAVE_BUTTON_WIDTH = 70
    save_ring_btn = Button(label="⬇ save", button_type="default", width=SAVE_BUTTON_WIDTH,
                            height=TOOLBAR_CONTROL_HEIGHT)
    save_zoom_btn = Button(label="⬇ save", button_type="default", width=SAVE_BUTTON_WIDTH,
                            height=TOOLBAR_CONTROL_HEIGHT)
    save_dotplot_btn = Button(label="⬇ save", button_type="default", width=SAVE_BUTTON_WIDTH,
                               height=TOOLBAR_CONTROL_HEIGHT)
    # Exports the cross blocks currently on screen (min block size filter
    # applied) as a links.tsv (docs/specs/hit_table.md section 5) -- lives on
    # the dotplot row since that's the panel showing every cross-genome block
    # at once, unlike the ring (also self-links) or the zoom panel (one
    # pivot/pair at a time). See SYN.exportBlocksTsv.
    export_blocks_btn = Button(label="⬇ blocks TSV", button_type="default", width=110,
                                height=TOOLBAR_CONTROL_HEIGHT)
    # format choice lives next to each save button rather than as a second
    # button per panel -- adds ~50px to a row instead of ~130px, which
    # matters on the dotplot row (already save button + order_toggle). JPEG
    # has no transparency channel -- SYN.exportFigureAsRaster fills white
    # first regardless of format, so this doesn't need special-casing there.
    EXPORT_FORMATS = ["SVG", "PNG", "JPEG"]
    ring_format_sel = Select(options=EXPORT_FORMATS, value="SVG", width=80,
                              height=TOOLBAR_CONTROL_HEIGHT)
    zoom_format_sel = Select(options=EXPORT_FORMATS, value="SVG", width=80,
                              height=TOOLBAR_CONTROL_HEIGHT)
    dotplot_format_sel = Select(options=EXPORT_FORMATS, value="SVG", width=80,
                                 height=TOOLBAR_CONTROL_HEIGHT)
    # off by default: natural (file/karyotype) order is what most users
    # recognize their chromosomes by, and reordering only pays off when the
    # two genomes are close enough that a near-1:1 correspondence exists to
    # reveal in the first place (see SYN.computeSimilarityOrder's docstring)
    order_toggle = Toggle(label="Order by similarity", active=False,
                           button_type="default", width=160, height=TOOLBAR_CONTROL_HEIGHT)
    # Global -- affects the ring AND the dotplot alike, but lives in the
    # ring's own row (next to show_gaps_toggle, itself global -- see its own
    # comment below) rather than the shared controls row above, since that's
    # where it was asked to sit.
    # On by default: matches today's only behavior (this plot script sorts
    # by size, by default -- see Dataset.__init__). When off, both panels
    # fall back to each
    # genome's own FASTA/natural sequence order. Layered under
    # order_toggle, not alongside it: while "order by similarity" is on,
    # this switch has no visible effect on the dotplot (similarity order
    # replaces natural order wholesale either way) -- it only governs what
    # "natural" itself means, for whichever panel is currently showing it.
    # See SYN.dpNaturalOrder (dotplot + the shared order concept) and
    # SYN.applyRingLayout (the ring, which has no similarity concept at all).
    # explicit width -- all four ring-row switches below now share one row
    # with the save button and format dropdown (see left_col), and none of
    # their unconstrained natural widths fit together under the ring's own
    # 620px. There's no width here that keeps every label on one line AND
    # the row under 620px -- Bokeh's Switch reserves real estate for the
    # toggle control itself before any label text, so even the shortest
    # label ("Show gaps") needs more room than this to stay on one line.
    # Trading that off deliberately: a label wrapping to two lines only
    # makes this one row a few px taller, not wider, so it doesn't reopen
    # the ring/zoom gap the row split (see left_col's git history) was
    # fixing in the first place -- width chosen for a comfortable margin
    # under 620px, not to chase single-line labels that don't fit regardless.
    RING_SWITCH_WIDTH = 95
    # tighter gap between each switch and its own label (Bokeh's own default
    # is 6px) and centered vertically against it -- matters especially now
    # that a wrapped two-line label (see RING_SWITCH_WIDTH above) is taller
    # than the toggle control itself, which otherwise sits pinned to the top
    # of that extra height instead of centered against the full label block.
    # :host is the styling entry point for every Bokeh widget's own shadow
    # DOM -- see Switch.stylesheets' docstring.
    RING_SWITCH_CSS = ':host{gap:2px;align-items:center;}'
    # same, plus a thin divider on this switch's leading edge, in the small
    # gap Bokeh already leaves between adjacent row items -- marks where one
    # text+toggle section ends and the next begins. Not applied to
    # self_links_toggle (the first of these four in the row -- see
    # left_col below): nothing of its own kind precedes it to divide from.
    # var(--divider-color) is Bokeh's own theme token for exactly this
    # (already used for the divider between this widget and its neighbors
    # in Bokeh's stock toolbars), not a hardcoded color of this file's own.
    RING_SWITCH_DIVIDER_CSS = (':host{gap:2px;align-items:center;'
                                'border-left:1px solid var(--divider-color);padding-left:6px;}')
    size_order_toggle = Switch(label="Order by size", active=True, width=RING_SWITCH_WIDTH,
                                stylesheets=[RING_SWITCH_DIVIDER_CSS])
    # Switch, not Toggle -- a checkbox-style on/off switch rather than a
    # pressable button, for these two specifically (order_toggle above stays
    # a button; only these ring-display toggles were asked to become
    # switches). off by default -- self-links can dominate/clutter the ring
    # (e.g. a heavily-homeologous polyploid genome), so a viewer opts in
    # rather than opting out. Lives with the ring's own save row (not the
    # palette/colors/min-block-size row below) because self-links are
    # drawn only on the ring -- SYN.buildDotplotSegmentsForLayout reads only
    # SYN.data.linksByQuery, which never contains homeolog records, so they
    # never appear on the zoom or dotplot panels either.
    self_links_toggle = Switch(label="Show self-links", active=False, width=RING_SWITCH_WIDTH,
                                stylesheets=[RING_SWITCH_CSS])
    # off by default -- hides cross-genome synteny ribbons, leaving only
    # self-links on screen if self_links_toggle is also on (see the combined
    # predicate in SYN.buildOverviewRibbons: is_homeolog records are governed
    # by self_links_toggle regardless of this switch; non-homeolog records
    # are governed by this switch regardless of self_links_toggle -- the two
    # compose independently, including the "both off" case, which is a valid
    # (if visually empty) combination, not specially prevented).
    hide_synteny_toggle = Switch(label="Hide synteny", active=False, width=RING_SWITCH_WIDTH,
                                  stylesheets=[RING_SWITCH_DIVIDER_CSS])
    # off by default -- an assembly gap marker is a diagnostic/QC detail
    # (see bin/rename_sequences.py's --out_gaps), not something every viewer
    # needs on screen by default, and a draft-quality assembly can have
    # thousands of
    # them at the default 100bp threshold. Lives in the ring's own row (the
    # user asked for it "below the ring plot") even though it also affects
    # the zoom panel and the dotplot (see SYN.applyGapVisibility/
    # SYN.applyDotplotGapVisibility/SYN.refreshDetail's gap handling) --
    # unlike self_links_toggle/hide_synteny_toggle, which are genuinely
    # ring-only, this one switch is intentionally the single on/off control
    # for every panel's gap markers, since "show gaps" is one concept
    # regardless of which panel is currently displaying them.
    show_gaps_toggle = Switch(label="Show gaps", active=False, width=RING_SWITCH_WIDTH,
                               stylesheets=[RING_SWITCH_DIVIDER_CSS])

    save_ring_btn.js_on_click(CustomJS(args=dict(fig=overview, fmt=ring_format_sel), code="""
        SYN.exportFigure(fig, 'ring', fmt.value);
    """))
    save_zoom_btn.js_on_click(CustomJS(args=dict(fig=detail_fig, fmt=zoom_format_sel), code="""
        SYN.exportFigure(fig, 'zoom', fmt.value);
    """))
    save_dotplot_btn.js_on_click(CustomJS(args=dict(fig=dotplot_fig, fmt=dotplot_format_sel), code="""
        SYN.exportFigure(fig, 'dotplot', fmt.value);
    """))
    export_blocks_btn.js_on_click(CustomJS(code="SYN.exportBlocksTsv();"))

    # Rebuilds every order-dependent dotplot source for whichever order is
    # now active -- SYN.applyDotplotOrder also updates SYN.state's offsets,
    # so the color/palette/min-block-size callbacks below (which redraw
    # the dotplot too, via SYN.applyDotplotSegmentsForCurrentLayout) pick up
    # the new order automatically without needing to know a reorder happened.
    order_toggle_callback = CustomJS(args=dict(
        color_spinner=color_spinner, min_block_spinner=min_block_spinner, dotplot_fig=dotplot_fig,
        dp_query_source=dp_q_src, dp_subject_source=dp_s_src, dp_grid_source=dp_grid_src,
        dp_q_label_source=dp_q_label_src, dp_s_label_source=dp_s_label_src,
        dp_segment_source=dp_seg_src, dp_gap_source=dp_gap_src,
    ), code="""
        // read by SYN.applyChainResult, so a fresh chain result knows
        // whether the dotplot needs a full similarity reorder or just a
        // segment redraw in place -- see that function's own comment
        SYN.state.orderBySimilarity = cb_obj.active;
        // SYN.dpOrderFor also applies the current min-length filter (see
        // min_seq_size_spinner below) on top of similarity/natural order, so
        // toggling this never silently drops that filter
        const order = SYN.dpOrderFor(cb_obj.active);
        SYN.applyDotplotOrder(order.queryOrder, order.subjectOrder, color_spinner.value, min_block_spinner.value, {
            query: dp_query_source, subject: dp_subject_source, grid: dp_grid_source,
            queryLabel: dp_q_label_source, subjectLabel: dp_s_label_source,
            segment: dp_segment_source, gap: dp_gap_source, fig: dotplot_fig,
        });
    """)
    order_toggle.js_on_click(order_toggle_callback)

    # size_order_toggle always reorders the ring (which has no similarity
    # concept to defer to); it only reorders the dotplot when order_toggle
    # is OFF, since while similarity order is active, "natural" isn't
    # currently visible on the dotplot at all -- see size_order_toggle's
    # own comment above for why this is layered rather than independent.
    size_order_toggle_callback = CustomJS(args=dict(
        query_source=q_src, subject_source=s_src, ribbon_source=r_src, label_source=label_src,
        color_spinner=color_spinner, min_block_spinner=min_block_spinner, dotplot_fig=dotplot_fig,
        self_links_toggle=self_links_toggle, hide_synteny_toggle=hide_synteny_toggle,
        order_toggle=order_toggle, gap_source=gap_src,
        dp_query_source=dp_q_src, dp_subject_source=dp_s_src, dp_grid_source=dp_grid_src,
        dp_q_label_source=dp_q_label_src, dp_s_label_source=dp_s_label_src,
        dp_segment_source=dp_seg_src, dp_gap_source=dp_gap_src,
    ), code="""
        SYN.state.orderBySize = cb_obj.active;
        // SYN.ringOrderFor/SYN.dpOrderFor apply the current min-length
        // filter (see min_seq_size_spinner below) on top of whichever
        // size/natural order this switch selects, so toggling it never
        // silently drops that filter.
        const ringOrder = SYN.ringOrderFor(cb_obj.active);
        SYN.applyRingLayout(ringOrder.queryOrder, ringOrder.subjectOrder, {
            querySource: query_source, subjectSource: subject_source,
            labelSource: label_source, ribbonSource: ribbon_source, gapSource: gap_source,
            minScore: min_block_spinner.value, k: color_spinner.value,
            showSelfLinks: self_links_toggle.active, hideSynteny: hide_synteny_toggle.active,
        });
        if (!order_toggle.active) {
            const dpOrder = SYN.dpOrderFor(false);
            SYN.applyDotplotOrder(dpOrder.queryOrder, dpOrder.subjectOrder, color_spinner.value, min_block_spinner.value, {
                query: dp_query_source, subject: dp_subject_source, grid: dp_grid_source,
                queryLabel: dp_q_label_source, subjectLabel: dp_s_label_source,
                segment: dp_segment_source, gap: dp_gap_source, fig: dotplot_fig,
            });
        }
    """)
    size_order_toggle.js_on_change('active', size_order_toggle_callback)

    # Global, like show_gaps_toggle above it -- filters both the ring and the
    # dotplot down to chromosomes at or above this length, composing with
    # whichever order/similarity toggle is currently active (via
    # SYN.ringOrderFor/SYN.dpOrderFor, the same helpers order_toggle/
    # size_order_toggle's own callbacks now go through). Any active pivot/
    # pair zoom is cleared on change, same as reset_btn -- the chromosome it
    # was zoomed into may no longer be visible at all.
    min_seq_size_callback = CustomJS(args=dict(
        order_toggle=order_toggle, size_order_toggle=size_order_toggle,
        color_spinner=color_spinner, min_block_spinner=min_block_spinner, dotplot_fig=dotplot_fig,
        query_source=q_src, subject_source=s_src, ribbon_source=r_src, label_source=label_src,
        gap_source=gap_src, self_links_toggle=self_links_toggle, hide_synteny_toggle=hide_synteny_toggle,
        dp_query_source=dp_q_src, dp_subject_source=dp_s_src, dp_grid_source=dp_grid_src,
        dp_q_label_source=dp_q_label_src, dp_s_label_source=dp_s_label_src,
        dp_segment_source=dp_seg_src, dp_gap_source=dp_gap_src,
        bar_source=detail_bar_src, detail_ribbon_source=detail_rib_src,
        detail_label_source=detail_label_src, detail_fig=detail_fig, detail_gap_source=detail_gap_src,
    ), code="""
        SYN.state.minSeqSize = cb_obj.value * 1e6;

        const ringOrder = SYN.ringOrderFor(size_order_toggle.active);
        SYN.applyRingLayout(ringOrder.queryOrder, ringOrder.subjectOrder, {
            querySource: query_source, subjectSource: subject_source,
            labelSource: label_source, ribbonSource: ribbon_source, gapSource: gap_source,
            minScore: min_block_spinner.value, k: color_spinner.value,
            showSelfLinks: self_links_toggle.active, hideSynteny: hide_synteny_toggle.active,
        });

        const dpOrder = SYN.dpOrderFor(order_toggle.active);
        SYN.applyDotplotOrder(dpOrder.queryOrder, dpOrder.subjectOrder, color_spinner.value, min_block_spinner.value, {
            query: dp_query_source, subject: dp_subject_source, grid: dp_grid_source,
            queryLabel: dp_q_label_source, subjectLabel: dp_s_label_source,
            segment: dp_segment_source, gap: dp_gap_source, fig: dotplot_fig,
        });

        // same as reset_btn's own callback -- a filtered-out chromosome may
        // be the current pivot/pair, and there's no cheap way to tell from
        // here, so always clear rather than risk a stale zoom panel
        SYN.state.mode = null;
        SYN.state.pivotSide = null;
        SYN.state.pivotName = null;
        SYN.state.pairTarget = null;
        SYN.state.pairSubject = null;
        query_source.selected.indices = [];
        subject_source.selected.indices = [];
        dp_query_source.selected.indices = [];
        dp_subject_source.selected.indices = [];
        ribbon_source.selected.indices = [];
        SYN.applyDetail(SYN.emptyDetail(), bar_source, detail_ribbon_source, detail_label_source,
                         detail_fig, detail_gap_source);
    """)
    min_seq_size_spinner.js_on_change('value', min_seq_size_callback)

    # Tapping a chromosome region -- a ring wedge or a dotplot ruler cell --
    # clears every OTHER clickable source's selection, so exactly one
    # chromosome (or pair) is ever "active" across all three panels at once.
    # other_sources is a list because there are four mutually-exclusive
    # .selected-driven click sources (ring x2, dotplot ruler x2), not two --
    # each needs to clear the other three (plus the ribbon source, a
    # different kind of selection -- see ribbon_tap_callback below). A
    # dotplot grid-square tap participates in the same exclusion too, but
    # through a plain Tap event handler instead of a source's own .selected
    # change -- see dotplot_fig's own Tap handler below for why (there is no
    # per-cell renderer/source left to select from any more). Two side
    # groups ('subject'/'query') rather than one callback per source because
    # cb_obj here is the Selection model that changed, not the
    # ColumnDataSource itself, so there's nothing reliable to branch on
    # inside a single shared callback.
    def make_tap_callback(this_source, other_sources, side):
        return CustomJS(args=dict(
            this_source=this_source, other_sources=other_sources, side=side,
            color_spinner=color_spinner, min_block_spinner=min_block_spinner, bar_source=detail_bar_src,
            ribbon_source=detail_rib_src, label_source=detail_label_src, detail_fig=detail_fig,
            detail_gap_source=detail_gap_src,
        ), code="""
            const idx = this_source.selected.indices;
            if (idx.length === 0) { return; }
            for (const other of other_sources) { other.selected.indices = []; }
            const name = this_source.data['name'][idx[idx.length - 1]];
            SYN.state.mode = 'pivot';
            SYN.state.pivotSide = side;
            SYN.state.pivotName = name;
            const d = SYN.buildDetailData(side, name, color_spinner.value, min_block_spinner.value);
            SYN.applyDetail(d, bar_source, ribbon_source, label_source, detail_fig, detail_gap_source);
        """)

    s_src.selected.js_on_change(
        'indices', make_tap_callback(s_src, [q_src, dp_s_src, dp_q_src, r_src], 'subject'))
    q_src.selected.js_on_change(
        'indices', make_tap_callback(q_src, [s_src, dp_s_src, dp_q_src, r_src], 'query'))
    dp_s_src.selected.js_on_change(
        'indices', make_tap_callback(dp_s_src, [s_src, q_src, dp_q_src, r_src], 'subject'))
    dp_q_src.selected.js_on_change(
        'indices', make_tap_callback(dp_q_src, [s_src, q_src, dp_s_src, r_src], 'query'))

    # Clicking a ribbon highlights it and dims the rest (SYN.applyRibbonHighlight)
    # instead of driving the zoom panel -- so this doesn't reuse
    # make_tap_callback above, which is specific to the pivot-zoom behavior.
    # Still participates in the same mutual-exclusion: a ribbon click clears
    # any active wedge/dotplot selection, and (via r_src's presence in THEIR
    # other_sources lists above) a wedge/dotplot click clears an active
    # ribbon highlight right back -- exactly one click source is ever
    # "active" across the whole ring+dotplot at a time.
    ribbon_tap_callback = CustomJS(args=dict(
        ribbon_source=r_src, other_sources=[s_src, q_src, dp_s_src, dp_q_src],
    ), code="""
        const idx = ribbon_source.selected.indices;
        if (idx.length > 0) {
            for (const other of other_sources) { other.selected.indices = []; }
        }
        SYN.applyRibbonHighlight(ribbon_source);
    """)
    r_src.selected.js_on_change('indices', ribbon_tap_callback)

    # Tapping a dotplot grid square zooms the detail panel into exactly that
    # (target chromosome, reference chromosome) pair -- unlike the pivot
    # callback above, this always shows exactly two bars, and handles zero
    # links as a normal result (see SYN.buildPairDetailData) rather than
    # refusing the click, since "no synteny here" is itself the answer for
    # an empty cell. A plain Tap event on the whole figure, not a renderer
    # selection -- there is no invisible per-cell renderer any more (see
    # module docstring: that was 48,841 hidden polygons for axolotl alone),
    # so SYN.dpCellAt hit-tests the tapped (x, y) by binary search instead;
    # it returns null for the ruler strips and the dead margin beyond either
    # axis, which keep their own TapTool-driven selection above, so this is
    # a no-op there rather than double-handling the same click.
    dotplot_pair_tap_callback = CustomJS(args=dict(
        other_sources=[s_src, q_src, dp_s_src, dp_q_src, r_src],
        color_spinner=color_spinner, min_block_spinner=min_block_spinner, bar_source=detail_bar_src,
        ribbon_source=detail_rib_src, label_source=detail_label_src, detail_fig=detail_fig,
        detail_gap_source=detail_gap_src,
    ), code="""
        const hit = SYN.dpCellAt(cb_obj.x, cb_obj.y);
        if (!hit) { return; }
        for (const other of other_sources) { other.selected.indices = []; }
        SYN.state.mode = 'pair';
        SYN.state.pairTarget = hit.targetName;
        SYN.state.pairSubject = hit.subjectName;
        const d = SYN.buildPairDetailData(hit.targetName, hit.subjectName, color_spinner.value, min_block_spinner.value);
        SYN.applyDetail(d, bar_source, ribbon_source, label_source, detail_fig, detail_gap_source);
    """)
    dotplot_fig.js_on_event(Tap, dotplot_pair_tap_callback)

    color_callback = CustomJS(args=dict(
        subject_source=s_src, dp_subject_source=dp_s_src,
        overview_ribbon_source=r_src, dotplot_segment_source=dp_seg_src, min_block_spinner=min_block_spinner,
        self_links_toggle=self_links_toggle, hide_synteny_toggle=hide_synteny_toggle,
        bar_source=detail_bar_src, ribbon_source=detail_rib_src,
        label_source=detail_label_src, detail_fig=detail_fig, detail_gap_source=detail_gap_src,
    ), code="""
        const k = cb_obj.value;
        SYN.recolorSubjectWedges(k, subject_source);
        SYN.recolorSubjectWedges(k, dp_subject_source);
        SYN.applyOverviewRibbons(min_block_spinner.value, k, self_links_toggle.active,
                                  hide_synteny_toggle.active, overview_ribbon_source);
        SYN.applyDotplotSegmentsForCurrentLayout(min_block_spinner.value, k, dotplot_segment_source);
        SYN.refreshDetail(k, min_block_spinner.value, bar_source, ribbon_source, label_source, detail_fig,
                           detail_gap_source);
    """)
    color_spinner.js_on_change('value', color_callback)

    min_block_callback = CustomJS(args=dict(
        overview_ribbon_source=r_src, dotplot_segment_source=dp_seg_src, color_spinner=color_spinner,
        self_links_toggle=self_links_toggle, hide_synteny_toggle=hide_synteny_toggle,
        bar_source=detail_bar_src, ribbon_source=detail_rib_src,
        label_source=detail_label_src, detail_fig=detail_fig, detail_gap_source=detail_gap_src,
    ), code="""
        const minScore = cb_obj.value;
        const k = color_spinner.value;
        SYN.applyOverviewRibbons(minScore, k, self_links_toggle.active,
                                  hide_synteny_toggle.active, overview_ribbon_source);
        SYN.applyDotplotSegmentsForCurrentLayout(minScore, k, dotplot_segment_source);
        SYN.refreshDetail(k, minScore, bar_source, ribbon_source, label_source, detail_fig, detail_gap_source);
        SYN.updateChainStatus();
    """)
    min_block_spinner.js_on_change('value', min_block_callback)

    # Min identity, max gap, and hit rank each need a fresh client-side
    # re-chain (SYN.requestChain) rather than a re-filter/re-render of
    # already-computed data -- unlike every control above. SYN.ui (set once
    # by SYN.init, before any of these can fire) already holds every widget
    # SYN.currentChainParams/SYN.applyChainResult need, so there's nothing
    # to pass as CustomJS args here.
    chain_params_callback = CustomJS(code="SYN.requestChain(SYN.currentChainParams());")
    min_identity_spinner.js_on_change('value', chain_params_callback)
    max_gap_spinner.js_on_change('value', chain_params_callback)
    hit_rank_select.js_on_change('value', chain_params_callback)

    # Switch fires via js_on_change('active', ...), not js_on_click (which is
    # Button/Toggle-specific) -- cb_obj is still the Switch itself either way,
    # so cb_obj.active works the same as it did when this was a Toggle.
    self_links_toggle_callback = CustomJS(args=dict(
        overview_ribbon_source=r_src, color_spinner=color_spinner, min_block_spinner=min_block_spinner,
        hide_synteny_toggle=hide_synteny_toggle,
    ), code="""
        SYN.applyOverviewRibbons(min_block_spinner.value, color_spinner.value, cb_obj.active,
                                  hide_synteny_toggle.active, overview_ribbon_source);
    """)
    self_links_toggle.js_on_change('active', self_links_toggle_callback)

    hide_synteny_toggle_callback = CustomJS(args=dict(
        overview_ribbon_source=r_src, color_spinner=color_spinner, min_block_spinner=min_block_spinner,
        self_links_toggle=self_links_toggle,
    ), code="""
        SYN.applyOverviewRibbons(min_block_spinner.value, color_spinner.value, self_links_toggle.active,
                                  cb_obj.active, overview_ribbon_source);
    """)
    hide_synteny_toggle.js_on_change('active', hide_synteny_toggle_callback)

    # Single on/off switch for every panel's gap markers at once (see
    # show_gaps_toggle's own comment above for why it's not split per
    # panel): flips SYN.state.showGaps, then rebuilds the ring's gap source
    # from the current-order master list (SYN.data.gapRecords), the
    # dotplot's gap lines from the current dotplot offsets, and re-renders
    # whatever's currently in the zoom panel (a pivot, a pair, or nothing --
    # SYN.refreshDetail is a no-op in the last case) so its gap ticks
    # (already computed into every SYN.buildDetailData/buildPairDetailData
    # result regardless of this switch -- see SYN.applyDetail) actually
    # show or hide too.
    show_gaps_toggle_callback = CustomJS(args=dict(
        gap_source=gap_src, dp_gap_source=dp_gap_src,
        color_spinner=color_spinner, min_block_spinner=min_block_spinner,
        bar_source=detail_bar_src, ribbon_source=detail_rib_src,
        label_source=detail_label_src, detail_fig=detail_fig, detail_gap_source=detail_gap_src,
    ), code="""
        SYN.state.showGaps = cb_obj.active;
        SYN.applyGapVisibility(gap_source);
        SYN.applyDotplotGapVisibility(dp_gap_source);
        SYN.refreshDetail(color_spinner.value, min_block_spinner.value, bar_source, ribbon_source,
                           label_source, detail_fig, detail_gap_source);
    """)
    show_gaps_toggle.js_on_change('active', show_gaps_toggle_callback)

    # Swaps the active color array itself (SYN.data.palette -- every recolor
    # function reads from it by reference, so nothing else needs to change)
    # rather than the number of colors cycled through, which is what
    # color_spinner controls -- the two are independent and compose (e.g.
    # "Pastel" at a count of 4).
    palette_callback = CustomJS(args=dict(
        subject_source=s_src, dp_subject_source=dp_s_src,
        overview_ribbon_source=r_src, dotplot_segment_source=dp_seg_src,
        color_spinner=color_spinner, min_block_spinner=min_block_spinner,
        self_links_toggle=self_links_toggle, hide_synteny_toggle=hide_synteny_toggle,
        bar_source=detail_bar_src, ribbon_source=detail_rib_src,
        label_source=detail_label_src, detail_fig=detail_fig, detail_gap_source=detail_gap_src,
    ), code="""
        SYN.data.palette = SYN.data.palettes[cb_obj.value];
        const k = color_spinner.value;
        SYN.recolorSubjectWedges(k, subject_source);
        SYN.recolorSubjectWedges(k, dp_subject_source);
        SYN.applyOverviewRibbons(min_block_spinner.value, k, self_links_toggle.active,
                                  hide_synteny_toggle.active, overview_ribbon_source);
        SYN.applyDotplotSegmentsForCurrentLayout(min_block_spinner.value, k, dotplot_segment_source);
        SYN.refreshDetail(k, min_block_spinner.value, bar_source, ribbon_source, label_source, detail_fig,
                           detail_gap_source);
    """)
    palette_select.js_on_change('value', palette_callback)

    reset_callback = CustomJS(args=dict(
        subject_source=s_src, query_source=q_src, dp_subject_source=dp_s_src, dp_query_source=dp_q_src,
        overview_ribbon_source=r_src,
        bar_source=detail_bar_src, ribbon_source=detail_rib_src,
        label_source=detail_label_src, detail_fig=detail_fig, detail_gap_source=detail_gap_src,
    ), code="""
        SYN.state.mode = null;
        SYN.state.pivotSide = null;
        SYN.state.pivotName = null;
        SYN.state.pairTarget = null;
        SYN.state.pairSubject = null;
        subject_source.selected.indices = [];
        query_source.selected.indices = [];
        dp_subject_source.selected.indices = [];
        dp_query_source.selected.indices = [];
        overview_ribbon_source.selected.indices = [];
        SYN.applyDetail(SYN.emptyDetail(), bar_source, ribbon_source, label_source, detail_fig, detail_gap_source);
    """)
    reset_btn.js_on_click(reset_callback)

    # not currently shown on the page -- left out of `layout` below, for now,
    # per explicit request to drop the usage instructions -- kept defined
    # here rather than deleted so it's a one-line change to bring back
    hint = Div(text="<p style='color:#666;font-size:13px'>Click any chromosome wedge on the "
                     "ring, or any chromosome band on the dotplot's axes, to zoom into its "
                     "links against the other genome -- both trigger the same zoom, and a "
                     "click on any of them clears the others' selection. You can also click a "
                     "single square in the dotplot grid to zoom straight into that one "
                     "chromosome pair, including squares with no alignments at all. Click a "
                     "ribbon on the ring instead to highlight just that one and dim the rest -- "
                     "click empty ring space, double-click the ring, or click a wedge/dotplot "
                     "cell to clear the highlight. Target "
                     "label/Reference label relabel the genomes everywhere a title or bar shows "
                     "them, so an exported panel can show the actual species/genome name "
                     "instead of \"target\"/\"reference\". Color palette switches the set of "
                     "colors reference chromosomes cycle through; Colors changes how many "
                     "discrete colors from it they cycle through. Min block size filters "
                     "out synteny blocks with fewer than that many supporting protein alignments "
                     "-- raise it to cut noise, lower it to see more (shorter, less certain) "
                     "blocks. Show self-links, next to the ring's save button, switches whether "
                     "a genome's own self-comparison links (e.g. homeologous chromosome pairs in "
                     "a polyploid genome) are drawn on the ring alongside the cross-genome "
                     "synteny -- off by default. Hide synteny, next to it, hides the cross-genome "
                     "links instead, leaving only self-links on screen if Show self-links is also "
                     "on -- off by default. The two switches are independent and compose. Show "
                     "gaps, next to those, marks assembly gaps (runs of N's in the input FASTA, "
                     "--min_asm_gap bp or longer) on the ring, the zoom panel, and (as thin dotted "
                     "lines) the dotplot alike -- off by default. Order "
                     "by size, next to those, sorts every chromosome, on the ring and both "
                     "dotplot axes, largest to smallest -- on by default; switch it off to see "
                     "them in their original FASTA order instead. Order chromosomes by similarity "
                     "reorders both dotplot axes (only -- the ring keeps whichever order the size "
                     "switch is set to) so shared synteny lines up into a diagonal instead -- most "
                     "useful when the two genomes are closely related with a roughly 1:1 "
                     "chromosome correspondence; off by default. Hover any "
                     "wedge, band, ribbon, or gap marker for details, and use each panel's save "
                     "button to "
                     "export it as SVG (a vector original -- open it in Inkscape, Illustrator, "
                     "or similar to edit it or convert it to PDF), PNG, or JPEG (both ready to "
                     "paste into a slide or document), whichever the dropdown next to that "
                     "button is set to.</p>")

    left_col = column(overview,
                       row(save_ring_btn, ring_format_sel, self_links_toggle, hide_synteny_toggle,
                           show_gaps_toggle, size_order_toggle))
    right_col = column(dotplot_fig,
                        row(save_dotplot_btn, dotplot_format_sel, export_blocks_btn, order_toggle))

    # static -- the tool's own name/tagline, not this run's target/reference
    # (that's the target/reference label inputs' job) -- unlike those, never
    # rewritten by label_callback
    header_title_div = Div(text="<h2>quick_synteny - fast and annotation-free synteny "
                                 "visualisation using spliced protein alignments</h2>")

    # No separate "target: <file>" / "reference: <file>" line here any more --
    # target_label_input/reference_label_input above default to exactly that
    # (query_subtitle/subject_subtitle, falling back to the role tag), so
    # showing the input file name a second time in its own Div would just be
    # the same information twice. That does mean the file name is no longer
    # visible once a viewer *retypes* a label -- accepted tradeoff, since the
    # box itself is still right there labeled "Target label"/"Reference label".

    # initial render matches target_label_input/reference_label_input's own
    # default (query_subtitle/subject_subtitle, falling back to the role
    # tag) -- reactive after that (see SYN.applyStats, wired below).
    # Empty/invisible if no --stats was given. width=STATS_PANEL_WIDTH makes
    # the box itself (format_stats_html/SYN.buildStatsHtml's markup, which is
    # `display:block` precisely so it stretches to fill this rather than
    # shrink-wrapping its content) as wide as the zoom panel's drawn frame
    # above it -- detail_fig's own width minus its toolbar strip (see
    # PLOT_TOOLBAR_WIDTH), not the full nominal figure width (which would run
    # this box's right edge out under the toolbar icons) and not the
    # narrower reset/save/format row under it either.
    # margin=(0,0,0,0) overrides Bokeh's own default 5px margin on every side
    # of a Div, which otherwise insets this box 5px from detail_fig's left
    # edge while pushing its right edge 5px past detail_fig's own -- exactly
    # the width but visibly not aligned with the panel it's meant to match.
    # The vertical breathing room that default margin used to provide is
    # already there regardless, from format_stats_html's own inline
    # `margin:4px 0 8px 0` on the box itself.
    stats_div = Div(text=format_stats_html(query_subtitle or query_name, subject_subtitle or subject_name,
                                            alignment_stats),
                     width=STATS_PANEL_WIDTH, margin=(0, 0, 0, 0))

    # placed under the zoom panel specifically (not a full-width bar under
    # all three panels) -- fills the otherwise-empty space under the
    # narrower middle column, since the ring/dotplot are taller
    mid_col = column(detail_fig, row(reset_btn, save_zoom_btn, zoom_format_sel), stats_div)

    label_callback = CustomJS(args=dict(
        target_label_input=target_label_input, reference_label_input=reference_label_input,
        overview=overview, dotplot_fig=dotplot_fig, stats_div=stats_div,
        q_src=q_src, dp_q_src=dp_q_src, s_src=s_src, dp_s_src=dp_s_src,
        color_spinner=color_spinner, min_block_spinner=min_block_spinner,
        bar_source=detail_bar_src, ribbon_source=detail_rib_src,
        label_source=detail_label_src, detail_fig=detail_fig, detail_gap_source=detail_gap_src,
    ), code="""
        const targetLabel = target_label_input.value.trim() || 'target';
        const referenceLabel = reference_label_input.value.trim() || 'reference';
        SYN.applyLabels(targetLabel, referenceLabel, {
            overview, dotplotFig: dotplot_fig,
            queryLikeSources: [q_src, dp_q_src], subjectLikeSources: [s_src, dp_s_src],
        });
        SYN.refreshDetail(color_spinner.value, min_block_spinner.value, bar_source, ribbon_source, label_source,
                           detail_fig, detail_gap_source);
        SYN.applyStats(targetLabel, referenceLabel, stats_div);
    """)
    target_label_input.js_on_change('value', label_callback)
    reference_label_input.js_on_change('value', label_callback)

    layout = column(
        header_title_div,
        row(reference_label_input, target_label_input, palette_select, color_spinner,
            min_identity_spinner, max_gap_spinner, hit_rank_select, min_block_spinner,
            min_seq_size_spinner, chain_status_div),
        # hint left out of the layout for now (not deleted -- still built
        # above, just not attached to anything file_html walks/serializes)
        # per explicit request to drop the usage instructions from the page
        # -- reinstate by adding `hint,` back here
        # spacing=20 -- row()'s default is 0, so without this the three
        # panels would sit flush against each other (or worse, apart by
        # whatever a child happens to overflow to, see left_col's own
        # comment above) rather than by a deliberate, equal gap
        row(left_col, mid_col, right_col, spacing=20),
    )

    # The page is built as an explicit Document (rather than handing `layout`
    # straight to file_html) so SYN.init can be wired to fire once, from the
    # browser's own "document ready" event, rather than needing Python to
    # pre-render any initial geometry (see module docstring). DocumentReady
    # fires after Bokeh has finished constructing every model client-side,
    # confirmed against both bokeh 3.9 and 3.10 (see this stream's plan) --
    # well before a viewer could touch a control, so SYN.init always wins
    # the race to set up SYN.state and every source's first render.
    doc = Document()
    doc.add_root(layout)
    # every widget/source SYN.init stashes into SYN.ui (see that function's
    # own comment) -- new code (SYN.requestChain, SYN.applyChainResult,
    # SYN.exportBlocksTsv) reads/writes through that registry instead of
    # each getting its own CustomJS args=dict(...), unlike the pre-existing
    # callbacks above, which are left as they were (see the plan's own
    # instruction not to refactor beyond what's needed).
    doc.js_on_event(DocumentReady, CustomJS(args=dict(
        q_src=q_src, s_src=s_src, r_src=r_src, label_src=label_src, gap_src=gap_src,
        dp_q_src=dp_q_src, dp_s_src=dp_s_src, dp_grid_src=dp_grid_src,
        dp_q_label_src=dp_q_label_src, dp_s_label_src=dp_s_label_src,
        dp_seg_src=dp_seg_src, dp_gap_src=dp_gap_src, dotplot_fig=dotplot_fig,
        color_spinner=color_spinner, min_block_spinner=min_block_spinner,
        min_identity_spinner=min_identity_spinner, max_gap_spinner=max_gap_spinner,
        hit_rank_select=hit_rank_select, self_links_toggle=self_links_toggle,
        hide_synteny_toggle=hide_synteny_toggle, chain_status_div=chain_status_div,
        bar_source=detail_bar_src, ribbon_source=detail_rib_src, label_source=detail_label_src,
        detail_fig=detail_fig, detail_gap_source=detail_gap_src,
    ), code="""
        SYN.init({
            querySource: q_src, subjectSource: s_src, ribbonSource: r_src,
            labelSource: label_src, gapSource: gap_src,
            dpQuerySource: dp_q_src, dpSubjectSource: dp_s_src, dpGridSource: dp_grid_src,
            dpQueryLabelSource: dp_q_label_src, dpSubjectLabelSource: dp_s_label_src,
            dpSegmentSource: dp_seg_src, dpGapSource: dp_gap_src, dotplotFig: dotplot_fig,
            colorSpinner: color_spinner, minBlockSpinner: min_block_spinner,
            minIdentitySpinner: min_identity_spinner, maxGapSpinner: max_gap_spinner,
            hitRankSelect: hit_rank_select, selfLinksToggle: self_links_toggle,
            hideSyntenyToggle: hide_synteny_toggle, chainStatusDiv: chain_status_div,
            detailBarSource: bar_source, detailRibbonSource: ribbon_source,
            detailLabelSource: label_source, detailFig: detail_fig, detailGapSource: detail_gap_source,
        });
    """))

    page_html = file_html(doc, CDN, title=f"{query_name} vs {subject_name} -- interactive synteny")

    syn_data = {
        'palette': PALETTES[DEFAULT_PALETTE_NAME],
        'palettes': PALETTES,
        'targetGrey': TARGET_GREY,
        'targetName': query_name,
        'referenceName': subject_name,
        # the *editable* label's own default (query_subtitle/subject_subtitle
        # -- the actual input file name/accession -- falling back to
        # "target"/"reference" when none was given, same fallback
        # target_label_input/reference_label_input themselves use, NOT
        # query_name/subject_name above -- those are "target"/"comparison",
        # the internal pipeline role tags (see main.nf), not this page's own
        # user-facing terms
        'targetLabelDefault': query_subtitle or query_name,
        'referenceLabelDefault': subject_subtitle or "reference",
        'queryNames': query_names,
        'subjectNames': subject_names,
        # natural (FASTA) order -- what SYN.dpNaturalOrder/SYN.buildRingLayout
        # fall back to when the "Order by size" switch is off. queryNames/
        # subjectNames above are already size order (see Dataset.__init__).
        'queryNamesNatural': [n for n, _ in ds.query_chroms_natural],
        'subjectNamesNatural': [n for n, _ in ds.subject_chroms_natural],
        'querySizes': ds.query_sizes,
        'subjectSizes': ds.subject_sizes,
        'subjectColorIndex': subject_index,
        # query (target) wedges are always flat TARGET_GREY, never palette-
        # indexed -- this index exists only so the target's own homeolog
        # ribbons can still be colored per-chromosome, same as
        # SYN.buildRingLayout's targetHomeologLinks loop
        'queryColorIndex': query_index,
        # Placeholders, replaced the moment the first client-side chain
        # result comes back (SYN.applyChainResult, called from
        # SYN.startChainer -- see SYN.init) -- Python computes no synteny at
        # all any more (see module docstring). Left as real (empty, not
        # missing) keys so SYN.init's own first ring/dotplot render, which
        # runs before that first chain result exists, has something safe to
        # read: zero ribbons/segments, not a crash.
        'linksByReference': {}, 'linksByQuery': {},
        'targetHomeologLinks': [], 'referenceHomeologLinks': [],
        'maxScore': 1,
        # per-chromosome assembly-gap lookup, one genome's own gaps each --
        # separate dicts (not a single by-name lookup) since target and
        # reference chromosomes can share names, same reasoning as
        # queryColorIndex/subjectColorIndex above. Read by SYN.buildGapRecords
        # (ring wedges), SYN.buildDotplotGapLines (dotplot lines), and
        # SYN.gapsForChrom (zoom panel bars) -- one source of raw gap data,
        # three different geometries built from it.
        'targetGapsByChrom': group_gaps_by_chrom(ds.query_gaps),
        'referenceGapsByChrom': group_gaps_by_chrom(ds.subject_gaps),
        # ring geometry constants -- embedded (not duplicated as separate JS
        # literals) so SYN.buildRingLayout can never numerically drift from
        # the fixed figure extent build_page() already committed to at
        # Python-render time (see OUTER_R/INNER_R/LINK_R/GROUP_GAP/CHROM_GAP
        # module constants above)
        'outerR': OUTER_R, 'innerR': INNER_R, 'linkR': LINK_R,
        'groupGap': GROUP_GAP, 'chromGap': CHROM_GAP, 'minGapAngle': MIN_GAP_ANGLE,
        'detailFigWidth': DETAIL_FIG_WIDTH,
        'statsPanelWidth': STATS_PANEL_WIDTH,
        'dpRulerFrac': DP_RULER_FRAC,
        # None (-> JSON null, falsy in JS) if --stats wasn't given -- every
        # reader of this (SYN.buildStatsHtml) already treats that as "no
        # panel to show"
        'alignmentStats': alignment_stats,
        # --min_identity/--min_block as given on the CLI, or None (-> JSON
        # null) meaning "use SYNCHAIN.autoParams" -- read once, by SYN.init,
        # the moment the embedded hit tables are decoded (see that
        # function's own comment). --max_gap always has a concrete value
        # (nextflow.config itself defaults it to 25), so max_gap_spinner's
        # own Python-side value= is already correct and this is only here
        # for SYN.currentChainParams' sake, not any JS-side "is it auto?"
        # branch the other two need.
        'initialParams': {'minPositive': min_identity, 'maxGap': max_gap, 'minBlock': min_block},
    }
    # bin/chain.js ahead of SHARED_JS (see module docstring): the one
    # implementation of the chaining algorithm, run in a Web Worker (or, if
    # that fails, synchronously on the main thread -- see SYN.startChainer).
    # A second, inert copy goes in its own <script type="text/plain"> tag
    # (BEFORE this one, so this stays the LAST <script> in the page --
    # tests/js/page_geometry.mjs and tests/js/page_chain_parity.mjs both
    # assume that) purely so SYN.startChainer can read the raw source text
    # back out at runtime to build the Worker's own Blob -- a Worker can't
    # share the main thread's already-parsed <script>, and this avoids a
    # runtime fetch()/XHR of bin/chain.js, which wouldn't work at all from a
    # file:// page (see the plan's own note on this).
    synchain_src_tag = f'<script type="text/plain" id="synchain-src">{chain_js_src}</script>\n'
    injected = (synchain_src_tag + "<script>\n" + chain_js_src + "\n" + SHARED_JS
                + "\nSYN.data = " + json.dumps(syn_data) + ";\nSYN.hitsPayload = "
                + json.dumps(hits_payload) + ";\n</script>\n</body>")
    return page_html.replace("</body>", injected, 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--query_name', required=True)
    parser.add_argument('--subject_name', required=True)
    parser.add_argument('--query_chrom_sizes', required=True)
    parser.add_argument('--subject_chrom_sizes', required=True)
    parser.add_argument('--target_hits', required=True,
                         help="bin/extract_hits.py output for the target genome "
                              "(docs/specs/hit_table.md section 1) -- embedded as SYN.hitsPayload "
                              "(section 6) and chained client-side; Python never chains anything")
    parser.add_argument('--comparison_hits', required=True, help='ditto, for the reference genome')
    parser.add_argument('--min_identity', type=float, default=None,
                         help='initial Min identity (%%) control value, as a 0-1 fraction -- '
                              'default (unset): SYNCHAIN.autoParams picks it client-side from the '
                              'actual hit tables')
    parser.add_argument('--max_gap', type=int, default=25,
                         help='initial Max gap (genes) control value (docs/specs/hit_table.md '
                              'section 3)')
    parser.add_argument('--min_block', type=int, default=None,
                         help='initial Min block size control value -- default (unset): '
                              'SYNCHAIN.autoParams picks it client-side, same as --min_identity')
    parser.add_argument('--stats', default=None,
                         help='optional compute_alignment_stats.py JSON -- shown as a small '
                              'always-visible summary panel on the page (proteome size, each '
                              "genome's aligned-protein count/mean identity) rather than baked "
                              'into any one exportable figure')
    parser.add_argument('--query_subtitle', default=None,
                         help='optional subtitle under the query/target name in the page '
                              'header -- e.g. the input file name, so a viewer can tell '
                              'which physical genome "target"/"reference" refer to')
    parser.add_argument('--subject_subtitle', default=None,
                         help='ditto, under the subject/reference name')
    parser.add_argument('--target_gaps', default=None,
                         help="optional bin/rename_sequences.py --out_gaps output (chrom, start, "
                              "end TSV, no header) for the target genome -- drawn as gap markers "
                              'on the ring/zoom panel/dotplot when the page\'s "Show gaps" switch is on')
    parser.add_argument('--comparison_gaps', default=None,
                         help='ditto, for the reference genome')
    parser.add_argument('--out_prefix', required=True)
    args = parser.parse_args()

    ds = Dataset(args.query_chrom_sizes, args.subject_chrom_sizes, args.target_gaps, args.comparison_gaps)
    if not ds.query_chroms or not ds.subject_chroms:
        sys.exit("ERROR: no chromosomes to plot -- check --min_seq_size isn't "
                  "filtering out everything")

    alignment_stats = None
    if args.stats:
        with open(args.stats) as f:
            alignment_stats = json.load(f)

    hits_payload = build_hits_payload(args.target_hits, args.comparison_hits,
                                       [n for n, _ in ds.query_chroms_natural],
                                       [n for n, _ in ds.subject_chroms_natural])
    chain_js_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'chain.js')
    with open(chain_js_path) as f:
        chain_js_src = f.read()

    page_html = build_page(
        ds, args.query_name, args.subject_name, args.query_subtitle, args.subject_subtitle,
        alignment_stats, hits_payload, chain_js_src, args.min_identity, args.max_gap, args.min_block)

    out_path = f"{args.out_prefix}.interactive.html"
    with open(out_path, 'w') as f:
        f.write(page_html)

    same_reference = hits_payload['genomes']['reference'] == 'same_as_target'
    target_n = hits_payload['genomes']['target']['n']
    reference_n = target_n if same_reference else hits_payload['genomes']['reference']['n']
    gap_bits = []
    if args.target_gaps or args.comparison_gaps:
        gap_bits.append(f"{len(ds.query_gaps)} target / {len(ds.subject_gaps)} reference assembly gap(s)")
    print(f"[plot_synteny_interactive] wrote {out_path} ({target_n} target / {reference_n} reference "
          f"hit(s) embedded, chained client-side"
          + ((f"; {'; '.join(gap_bits)}") if gap_bits else '') + ")",
          file=sys.stderr)


if __name__ == '__main__':
    main()
