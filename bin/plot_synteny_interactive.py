#!/usr/bin/env python3
"""Render a two-genome synteny plot as a self-contained interactive HTML page
(Bokeh, CustomJS only, BokehJS inlined -- no server, works offline). This is
the pipeline's only plot output; each panel exports itself as SVG/PNG
(SYN.exportFigure).

Three panels: a Circos-style ring, a linear zoom panel, and a whole-genome
dotplot. Clicking a chromosome (ring wedge or dotplot ruler) zooms the middle
panel into it, with its partners on the other side packed alongside; clicking
a dotplot cell zooms into that one chromosome pair (the only way to select a
pair with zero links, see SYN.buildPairDetailData). Reference chromosomes are
always drawn on top.

Python computes no geometry and no synteny blocks: it ships each genome's raw
hit table (bin/chain.js's EMBEDDED PAYLOAD header, SYN.hitsPayload) plus
chromosome sizes and palette indices as embedded JSON, and inlines
bin/chain.js so the browser chains client-side (in a Web Worker when
available), re-chaining on every min identity / max gap / hit rank change
(SYN.requestChain). Min block size is a pure client-side filter: extraction
doesn't depend on it (bin/chain.js's CHAINING RULES header, rule 6), so the
chainer always runs at minBlock=3. SYN.init renders the initial page through
the same SHARED_JS functions every later control change reuses.

Needs Bokeh, so it runs on envs/bokeh.yml (see
modules/local/plot_synteny.nf).
"""
import argparse
import base64
import gzip
import json
import math
import os
import re
import sys

import numpy as np

from bokeh.document import Document
from bokeh.embed import file_html
from bokeh.events import DocumentReady, DoubleTap, Reset, Tap
from bokeh.layouts import column, row
from bokeh.models import (ColumnDataSource, HoverTool, TapTool, CustomJS, CustomJSTickFormatter,
                           Button, Div, InlineStyleSheet, Label, Range1d, Select, Spinner, Switch,
                           TextInput, Tooltip, UIElement)
from bokeh.models.dom import HTML
from bokeh.plotting import figure
from bokeh.resources import INLINE

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
# seaborn's "dark" set -- same hues again, deepened/desaturated
PALETTE_DARK = [
    '#001C7F', '#B1400D', '#12711C', '#8C0800', '#591E71',
    '#592F0D', '#A23582', '#3C3C3C', '#B8850A', '#006374',
]
# matplotlib's "tab10" -- distinct hues from the seaborn sets above
PALETTE_TABLEAU = [
    '#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD',
    '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF',
]
# every option needs exactly MAX_COLORS entries: SYN.paletteColor indexes
# with idx % k for any k up to MAX_COLORS
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
# the target genome's self-links (homeologs) -- darker than TARGET_GREY so
# they stay visible against the target's own grey wedges, and never a
# palette color, which would read as a reference chromosome's
TARGET_SELF_GREY = '#666666'

# Non-selection glyph alpha while something is selected (SYN.select). A
# single selected block dims nearly everything, so ribbons/segments can't go
# much lower without the panel going blank; wedges/ruler bands dim less so
# the ring and axes stay readable.
RIBBON_DIM_ALPHA = 0.2
WEDGE_DIM_ALPHA = 0.5

OUTER_R = 1.0
RING_WIDTH = 0.045
INNER_R = OUTER_R - RING_WIDTH
LINK_R = INNER_R
GROUP_GAP = math.radians(4)
CHROM_GAP = math.radians(0.5)
# a ring wedge narrower than this gets no name label (hover/click still
# work): at the label radius 2.5 degrees is ~12 px of arc, about one line of
# the 10 px label text, so neighbouring labels can't pile up. In data space,
# so zooming the ring doesn't bring hidden labels back.
RING_LABEL_MIN_ANGLE = math.radians(2.5)

# minimum drawn angular width of an assembly gap on the ring: at true width
# most gaps (a few hundred bp) would be sub-pixel and unhoverable. The
# tooltip still shows the real length. Small on purpose: a gap wedge spans
# the ideogram's full radial height, so it should read as a thin needle.
MIN_GAP_ANGLE = math.radians(0.2)

# the zoom panel's pixel width -- embedded so SYN.widenGapBp's px-to-bp
# conversion matches the actual figure
DETAIL_FIG_WIDTH = 446

# Bokeh's right-hand toolbar strip, which sits inside the figure's width
PLOT_TOOLBAR_WIDTH = 30

# detail_fig's min_border_left, set explicitly so the stats box can line up
# with the drawn frame
DETAIL_FRAME_LEFT = 5

# the detail panel's two rows, in data units: reference bars span
# DETAIL_TOP_Y..+DETAIL_BAR_H, target bars DETAIL_BOT_Y..+DETAIL_BAR_H,
# ribbons in between (embedded as SYN.data.detailRows)
DETAIL_BAR_H = 0.32
DETAIL_TOP_Y = 1.2
DETAIL_BOT_Y = 0.0

# the stats box spans the zoom panel's drawn frame, not its nominal width
STATS_PANEL_WIDTH = DETAIL_FIG_WIDTH - PLOT_TOOLBAR_WIDTH - DETAIL_FRAME_LEFT

# shared height for every titleless toolbar control: left unset, buttons and
# Selects auto-size to their own label/font metrics and don't line up
TOOLBAR_CONTROL_HEIGHT = 32

# the header row's titled controls
TOP_CONTROL_WIDTH = 150

# The desktop app's light tokens (gui/renderer/style.css). Figures stay
# white so exports are unchanged; no dark mode for the same reason. Bokeh
# draws each widget in its own shadow root, which page CSS can't reach, and
# redefines its CSS variables on each :host -- so THEME_CSS goes on every
# element (see build_page()), after Bokeh's own.
UI_FONT = "system-ui, -apple-system, 'Segoe UI', sans-serif"
THEME = {
    'ground': '#F4F5F1', 'card': '#FFFFFF', 'line': '#E2E4DC', 'field_line': '#C9CCC0',
    'fg': '#1C1E19', 'fg2': '#5A5F55', 'fg3': '#6C7166',
    'accent': '#95B540', 'accent_ink': '#46600F', 'accent_tint': '#EEF4DE', 'accent_tint_fg': '#3A500C',
}
THEME_CSS = f""":host {{
    --bokeh-base-font: {UI_FONT}; --base-font: {UI_FONT};
    --color: {THEME['fg']}; --background-color: {THEME['card']}; --hover-color: {THEME['ground']};
    --border-color: {THEME['field_line']}; --divider-color: {THEME['line']}; --border-radius: 8px;
    --highlight-color: {THEME['accent']}; --active-border: {THEME['accent']};
    --active-bg: #C5DB8C; --active-fg: {THEME['accent_ink']};
    --inactive-bg: {THEME['line']}; --inactive-fg: #9EA393;
    --input-focus-border-color: {THEME['accent']}; --input-focus-halo-color: rgba(149, 181, 64, 0.35);
    --icon-color: {THEME['fg3']}; --shortcut-color: {THEME['fg3']}; --placeholder-color: #8E9386;
    --surface-background-color: #F1F2EC;
    --default: {THEME['card']}; --default-border: {THEME['line']};
    --default-hover: {THEME['ground']}; --default-hover-border: {THEME['field_line']};
    --default-active: {THEME['accent_tint']}; --default-active-border: {THEME['accent']};
}}
label {{ color: {THEME['fg2']}; }}"""
# a panel or the controls, as the app draws a card
CARD_CSS = (f":host{{background:{THEME['card']};border:1px solid {THEME['line']};border-radius:14px;"
            f"padding:8px 10px;box-shadow:0 1px 3px rgba(30,40,10,0.08);}}")
# the page around the cards; `html body` outranks the `html, body` reset
# file_html puts after it
PAGE_CSS = (f"<style>:root{{color-scheme:light;}}html body{{margin:0;padding:20px 16px 32px;"
            f"background:{THEME['ground']};color:{THEME['fg']};font-family:{UI_FONT};}}</style>")

# "Chromosome order" options -> (order by size, order by similarity), as
# SYN.data.orderModes. Similarity starts from size order (see
# SYN.computeSimilarityOrder), so chromosomes with no visible blocks tiebreak
# by size.
ORDER_MODES = {'Size': (True, False), 'File order': (False, False), 'Similarity': (True, True)}

# dotplot gap lines only -- lighter than the ring/zoom panel's black markers
GAP_LINE_COLOR = '#AAAAAA'

MAX_COLORS = 10
DEFAULT_COLORS = MAX_COLORS
# the "Min sequence length (Mb)" spinner's floor and starting value -- the
# pipeline's default --min_seq_size (10 kb)
MIN_SEQ_SIZE_FLOOR_MB = 0.01

# fraction of each dotplot axis given to the clickable chromosome-ruler
# strips that stand in for tick labels (Bokeh has no clickable ticks)
DP_RULER_FRAC = 0.035

ACCESSION_RE = re.compile(r'GC[AF]_\d+\.\d+')
FASTA_SUFFIX_RE = re.compile(r'(_genomic)?\.(fa|fna|fasta|fas)(\.b?gz)?$', re.IGNORECASE)


def short_species(species):
    """'Drosophila melanogaster' -> 'D. melanogaster': abbreviates the genus
    of a binomial (or longer) name, leaving anything that doesn't look like
    one (a single word, a lowercase first word) as it is."""
    words = species.split()
    if len(words) >= 2 and words[0][:1].isupper():
        return ' '.join([words[0][0] + '.'] + words[1:])
    return species


def source_label(source):
    """The shortest recognizable name for an input genome's source: its NCBI
    accession if it has one, else its file name minus the FASTA/compression
    extensions."""
    if not source:
        return None
    m = ACCESSION_RE.search(source)
    if m:
        return m.group(0)
    return FASTA_SUFFIX_RE.sub('', os.path.basename(source)) or source


def default_labels(stats, target_subtitle, reference_subtitle, target_fallback, reference_fallback):
    """The Target/Reference label inputs' starting values: short species
    names when compute_alignment_stats.py knows them ('D. melanogaster'),
    plus each accession when both are the same species; else the
    accession/file name, else the role word."""
    stats = stats or {}
    sides = [(stats.get('query_species'), target_subtitle or stats.get('query_source'), target_fallback),
             (stats.get('subject_species'), reference_subtitle or stats.get('subject_source'), reference_fallback)]
    same_species = bool(sides[0][0]) and sides[0][0] == sides[1][0]
    labels = []
    for species, source, fallback in sides:
        src = source_label(source)
        if species:
            label = short_species(species)
            labels.append(f"{label} {src}" if same_species and src else label)
        else:
            labels.append(src or fallback)
    return labels[0], labels[1]


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


def read_original_names(path, sizes):
    """bin/rename_sequences.py's --out_lookup output (original_id, new_id,
    length, type; one header line) -> {new_id: original_id}, for the
    chromosomes in `sizes` whose name actually changed. {} for path=None."""
    if not path:
        return {}
    names = {}
    with open(path) as f:
        next(f, None)  # header
        for line in f:
            if not line.strip():
                continue
            original_id, new_id = line.rstrip('\n').split('\t')[:2]
            if new_id in sizes and new_id != original_id:
                names[new_id] = original_id
    return names


def read_gaps(path):
    """bin/rename_sequences.py's --out_gaps output: chrom, start, end -- no
    header, 0-based half-open. Returns [] for path=None (no --target_gaps/
    --reference_gaps given)."""
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
    """Parses one gzipped bin/extract_hits.py output into per-column lists,
    in TSV row order (sorted by chrom, start, end, protein, rank).
    `identity`/`score` aren't part of chain.js's HitTable and are dropped."""
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
    """One GenomeColumns column (bin/chain.js's EMBEDDED PAYLOAD header):
    little-endian typed array, gzip'd, base64'd."""
    np_dtype = {'u8': '<u1', 'u16': '<u2', 'u32': '<u4', 'f64': '<f8'}[dtype]
    data = np.asarray(values, dtype=np_dtype).tobytes()
    return {'dtype': dtype, 'data': base64.b64encode(gzip.compress(data, mtime=0)).decode('ascii')}


def widened_dtype(values, narrow, wide, limit):
    """`wide` if any value reaches `limit` (startDelta/length go f64 past
    2**32, per bin/chain.js's EMBEDDED PAYLOAD header)."""
    return wide if values and max(values) >= limit else narrow


def build_genome_columns(hits, chrom_names, prot_index):
    """One GenomeColumns dict (bin/chain.js's EMBEDDED PAYLOAD header) for a
    parsed hit table. chrom_names is the genome's full, unfiltered
    chrom.sizes-order list; bin/extract_hits.py already restricted the hits
    to those sequences."""
    chrom_index = {name: i for i, name in enumerate(chrom_names)}
    n = len(hits['chrom'])
    chrom_idx = [chrom_index[c] for c in hits['chrom']]

    # start minus the previous row's start on the same chrom (first row per
    # chrom relative to 0); rows are already grouped by chrom
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


def build_hits_payload(target_hits_path, reference_hits_path, target_chrom_names, reference_chrom_names):
    """SYN.hitsPayload (bin/chain.js's EMBEDDED PAYLOAD header). The
    reference is `'same_as_target'` only when both the hit tables and the
    chrom.sizes orderings are identical: decodePayload reuses the target's
    chrom-index column under the reference's chromNames, which is only
    correct if chromosome N means the same thing on both sides."""
    with open(target_hits_path, 'rb') as f:
        target_bytes = f.read()
    same_genome = target_chrom_names == reference_chrom_names
    if same_genome:
        with open(reference_hits_path, 'rb') as f:
            same_genome = f.read() == target_bytes

    target_hits = read_hits_tsv(target_bytes)
    reference_hits = None if same_genome else read_hits_tsv(open(reference_hits_path, 'rb').read())

    # sorted union of both genomes' accessions (bin/chain.js's HIT TABLE
    # header); ASCII sorts the same in Python and JS's default Array.sort()
    all_proteins = set(target_hits['protein']) | (set(reference_hits['protein']) if reference_hits else set())
    proteins = sorted(all_proteins)
    prot_index = {p: i for i, p in enumerate(proteins)}

    target_cols = build_genome_columns(target_hits, target_chrom_names, prot_index)
    reference_cols = ('same_as_target' if same_genome
                       else build_genome_columns(reference_hits, reference_chrom_names, prot_index))

    proteins_gz = gzip.compress('\n'.join(proteins).encode('utf-8'), mtime=0)
    return {
        # bin/chain.js's VERSION: decodePayload refuses any other
        'version': 2,
        'proteins': base64.b64encode(proteins_gz).decode('ascii'),
        'genomes': {'target': target_cols, 'reference': reference_cols},
    }


class Dataset:
    def __init__(self, query_chrom_sizes, subject_chrom_sizes, query_gaps_path=None, subject_gaps_path=None,
                 query_lookup_path=None, subject_lookup_path=None):
        # chrom.sizes is in natural (FASTA) order. query_chroms/subject_chroms
        # are size order, matching order_select's "Size" default; the natural
        # lists are embedded for the "File order" mode and as the full,
        # unfiltered lists the hit tables' chrom indices refer to.
        self.query_chroms_natural = read_chrom_sizes(query_chrom_sizes)
        self.subject_chroms_natural = read_chrom_sizes(subject_chrom_sizes)
        self.query_chroms = sorted(self.query_chroms_natural, key=lambda c: -c[1])
        self.subject_chroms = sorted(self.subject_chroms_natural, key=lambda c: -c[1])
        self.query_sizes = dict(self.query_chroms)
        self.subject_sizes = dict(self.subject_chroms)

        # assembly gaps, restricted to chroms actually present
        self.query_gaps = [g for g in read_gaps(query_gaps_path) if g['chrom'] in self.query_sizes]
        self.subject_gaps = [g for g in read_gaps(subject_gaps_path) if g['chrom'] in self.subject_sizes]

        self.query_original_names = read_original_names(query_lookup_path, self.query_sizes)
        self.subject_original_names = read_original_names(subject_lookup_path, self.subject_sizes)


def group_gaps_by_chrom(gaps):
    """{chrom: [{start, end}, ...]}, embedded once and shared by all three
    panels (SYN.data.targetGapsByChrom/referenceGapsByChrom)."""
    by_chrom = {}
    for g in gaps:
        by_chrom.setdefault(g['chrom'], []).append({'start': g['start'], 'end': g['end']})
    return by_chrom


def build_overview_sources():
    """The ring's ColumnDataSources (query/subject wedges, ribbons, self
    ribbons, labels), with their columns but no rows -- SYN.init fills them
    via SYN.applyRingLayout on DocumentReady."""
    q_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[], name=[], size_label=[], group=[], display_name=[]))
    s_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[], name=[], size_label=[], group=[], palette_index=[],
                                  display_name=[]))
    r_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[], alpha=[], label=[], palette_index=[],
                                  q_chrom=[], s_chrom=[], key=[]))
    # self-links (homeologs) get their own renderer, styled apart from r_src
    r_self_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[], line_color=[], alpha=[], label=[],
                                       q_chrom=[], s_chrom=[], key=[], genome=[]))
    label_src = ColumnDataSource(dict(x=[], y=[], text=[]))
    return q_src, s_src, r_src, r_self_src, label_src


def build_dotplot_sources():
    """The dotplot's ColumnDataSources (query/subject ruler strips, block
    segments, gridlines, axis labels), with their columns but no rows --
    SYN.init fills them via SYN.applyDotplotOrder, which also sets the
    figure's x_range/y_range."""
    dp_q_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[], name=[], size_label=[], group=[], display_name=[]))
    dp_s_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[], name=[], size_label=[], group=[], palette_index=[],
                                  display_name=[]))
    dp_seg_src = ColumnDataSource(dict(xs=[], ys=[], line_color=[], alpha=[], label=[], palette_index=[],
                                       q_chrom=[], s_chrom=[], key=[]))
    dp_grid_src = ColumnDataSource(dict(xs=[], ys=[]))
    dp_q_label_src = ColumnDataSource(dict(x=[], y=[], text=[]))
    dp_s_label_src = ColumnDataSource(dict(x=[], y=[], text=[]))
    return dp_q_src, dp_s_src, dp_seg_src, dp_grid_src, dp_q_label_src, dp_s_label_src


# All client-side layout, chaining glue, and interaction (see module docstring).
SHARED_JS = r"""
window.SYN = window.SYN || {};
SYN.state = {
    // the one thing selected across all three panels (SYN.select), or null
    selection: null,
    // the min block size the viewer set, and the one SYN.applyChainResult
    // last fitted it to (a new longest block), until its callback has run
    minBlockWanted: null,
    minBlockFitted: null,
    // panel -> what the click being handled hit (SYN.onTapHit/SYN.resolveTap)
    pendingTap: {},
    // the dotplot's current per-chromosome axis offsets and orders (set only
    // by SYN.applyDotplotOrder); redraws that aren't reorders reuse them, and
    // SYN.dpCellAt binary-searches the orders to hit-test a tap
    dpQueryOffsets: null, dpSubjectOffsets: null,
    dpQueryOrder: null, dpSubjectOrder: null,
    // viewer-editable display names (SYN.applyLabels)
    targetLabel: null, referenceLabel: null,
    // the zoom panel's ranges change on every click, so its double-click
    // reset tracks the last programmatic range (SYN.applyDetail)
    detailRange: null,
    // the "Chromosome order" menu as two flags (ORDER_MODES): size vs FASTA
    // order (SYN.dpNaturalOrder), and similarity order on top, which makes
    // every re-chain and min-block-size change re-order. Defaults must match
    // order_select's ("Size"), as must the other Python/JS default pairs here.
    orderBySize: true,
    orderBySimilarity: false,
    // matches show_gaps_toggle's default
    showGaps: false,
    // the "Original names" switch (SYN.displayName) -- display only; plot
    // names stay the keys everywhere
    originalNames: false,
    // bp floor read by SYN.filterBySize; starts at min_seq_size_spinner's
    // floor (MIN_SEQ_SIZE_FLOOR_MB), so a default run hides nothing
    minSeqSize: 10000,
    // the dotplot's total span and ruler-strip thickness, which the min-length
    // filter can change; refreshed by SYN.applyDotplotOrder
    dpTotalX: null, dpTotalY: null, dpRx: null, dpRy: null,
};

// The label inputs are free text, interpolated into raw HTML (Div .text and
// HoverTool "{safe}" fields aren't escaped by Bokeh), so escape them first.
SYN.escapeHtml = function(s) {
    return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;').replace(/'/g, '&#39;');
};

// A block's look in every panel: its reference chromosome's palette colour,
// lightened by LINK_LIGHTEN, at an opacity rising with its anchor count from
// LINK_ALPHA_MIN to LINK_ALPHA_MAX (the largest block on the page, or in the
// detail panel)
SYN.LINK_LIGHTEN = 0.2;
SYN.LINK_ALPHA_MIN = 0.35;
SYN.LINK_ALPHA_MAX = 0.9;
SYN.linkAlpha = function(score, maxScore) {
    return SYN.LINK_ALPHA_MIN + (SYN.LINK_ALPHA_MAX - SYN.LINK_ALPHA_MIN) * Math.min(1, score / maxScore);
};
SYN.linkColor = function(hex) { return SYN.lighten(hex, SYN.LINK_LIGHTEN); };

SYN.lighten = function(hex, amount) {
    amount = amount === undefined ? 0.35 : amount;
    const r = parseInt(hex.slice(1, 3), 16), g = parseInt(hex.slice(3, 5), 16), b = parseInt(hex.slice(5, 7), 16);
    const blend = (c) => Math.round(c + (255 - c) * amount);
    const toHex = (c) => c.toString(16).padStart(2, '0');
    return '#' + toHex(blend(r)) + toHex(blend(g)) + toHex(blend(b));
};

SYN.darken = function(hex, amount) {
    amount = amount === undefined ? 0.35 : amount;
    const r = parseInt(hex.slice(1, 3), 16), g = parseInt(hex.slice(3, 5), 16), b = parseInt(hex.slice(5, 7), 16);
    const scale = (c) => Math.round(c * (1 - amount));
    const toHex = (c) => c.toString(16).padStart(2, '0');
    return '#' + toHex(scale(r)) + toHex(scale(g)) + toHex(scale(b));
};

// a block's identity across every source that draws it (ring ribbon,
// dotplot segment, detail ribbon) -- what a block selection stores
SYN.linkKey = function(l) {
    return `${l.q_chrom}:${l.q_start}-${l.q_end}|${l.s_chrom}:${l.s_start}-${l.s_end}|${l.orientation}`;
};

SYN.paletteColor = function(index, k) {
    return SYN.data.palette[index % k];
};

// black or white label text for a bar fill (some palettes include
// near-black swatches); simplified relative luminance, no gamma correction
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

// SYN.data.ribbons holds every ribbon (the chainer always runs at
// minBlock=3, the least strict floor bin/chain.js's rule 6 allows), so the
// min-block filter, recoloring and self-link visibility are all rebuilt
// together from that master list rather than by mutating what's displayed,
// which would let them undo each other.
//
// Returns {cross, self} for separate renderers, so self-links (is_homeolog: a
// protein hitting two loci in the same genome) can be drawn dashed over a
// fainter fill. A target self-link is dark grey rather than a palette color,
// which would read as a reference chromosome's.
SYN.buildOverviewRibbons = function(minScore, k, showSelfLinks, showSynteny) {
    const cross = SYN.data.ribbons.filter((r) => !r.is_homeolog && showSynteny && r.score >= minScore);
    const self = SYN.data.ribbons.filter((r) => r.is_homeolog && showSelfLinks && r.score >= minScore);
    const selfColor = (r) => r.genome === 'target' ? SYN.data.targetSelfGrey : SYN.paletteColor(r.palette_index, k);
    return {
        cross: {
            xs: cross.map((r) => r.xs),
            ys: cross.map((r) => r.ys),
            fill_color: cross.map((r) => SYN.linkColor(SYN.paletteColor(r.palette_index, k))),
            alpha: cross.map((r) => r.alpha),
            label: cross.map((r) => r.label),
            palette_index: cross.map((r) => r.palette_index),
            q_chrom: cross.map((r) => r.q_chrom), s_chrom: cross.map((r) => r.s_chrom), key: cross.map((r) => r.key),
        },
        self: {
            xs: self.map((r) => r.xs),
            ys: self.map((r) => r.ys),
            fill_color: self.map((r) => SYN.linkColor(selfColor(r))),
            line_color: self.map((r) => SYN.darken(selfColor(r), 0.2)),
            alpha: self.map((r) => r.alpha * SYN.SELF_LINK_ALPHA_SCALE),
            label: self.map((r) => r.label),
            q_chrom: self.map((r) => r.q_chrom), s_chrom: self.map((r) => r.s_chrom), key: self.map((r) => r.key),
            genome: self.map((r) => r.genome),
        },
    };
};

// self-link fill alpha relative to a cross ribbon of the same score, so a
// polyploid's many self-links don't drown out the synteny
SYN.SELF_LINK_ALPHA_SCALE = 0.6;

// ribbonSource is the cross-genome one; the self-link source always comes
// from SYN.ui (set by SYN.init before anything can call this)
SYN.applyOverviewRibbons = function(minScore, k, showSelfLinks, showSynteny, ribbonSource) {
    const built = SYN.buildOverviewRibbons(minScore, k, showSelfLinks, showSynteny);
    ribbonSource.data = built.cross;
    SYN.ui.selfRibbonSource.data = built.self;
    SYN.ui.selfRibbonSource.change.emit();
    ribbonSource.change.emit();
    // rebuilt rows no longer match the old selection indices
    SYN.reapplySelection();
};

// ---- ring geometry, up through SYN.buildRingLayout ----

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

// Like SYN.arcPoints, but walks a0 -> a1 as given, never reordered:
// SYN.ribbonPolygonJS needs to know which endpoint is which.
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

// An inverted ('-') block joins target start to reference end: swapping the
// reference span's ends twists the ribbon, marking the inversion.
SYN.orientEnds = function(s1, s2, link) {
    return link.orientation === '-' ? [s2, s1] : [s1, s2];
};

SYN.ribbonPolygonJS = function(aq1, aq2, as1, as2, radius, n) {
    n = n || 24;
    // query arc, Bezier through the centre to the subject arc, subject arc,
    // Bezier back
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

// min/max, since a reference chromosome's angle decreases with bp
SYN.widenGapAngle = function(a0, a1) {
    let lo = Math.min(a0, a1), hi = Math.max(a0, a1);
    if (hi - lo < SYN.data.minGapAngle) {
        const mid = (lo + hi) / 2;
        lo = mid - SYN.data.minGapAngle / 2;
        hi = mid + SYN.data.minGapAngle / 2;
    }
    return [lo, hi];
};

// Zoom-panel counterpart of SYN.widenGapAngle: floors a gap tick to
// GAP_MIN_PX. refSpan is the caller's estimate of the panel's data width;
// 1.08 is the detail x_range's span (-0.03 to 1.05 of xMax).
SYN.GAP_MIN_PX = 1.5;
SYN.widenGapBp = function(start, end, refSpan) {
    const minBp = (refSpan * 1.08 / SYN.data.detailFigWidth) * SYN.GAP_MIN_PX;
    if (end - start < minBp) {
        const mid = (start + end) / 2;
        return [mid - minBp / 2, mid + minBp / 2];
    }
    return [start, end];
};

// Ribbon master list (SYN.data.ribbons) for the given ring offsets -- split
// out of SYN.buildRingLayout so a re-chain can rebuild just the ribbons (on a
// gappy draft assembly the gap marks are the costliest part of a layout).
SYN.buildRibbonRecords = function(queryOffsets, subjectOffsets) {
    // each group sorted by score ascending, so higher-scoring ribbons are
    // drawn on top
    const ribbons = [];
    const allCross = [].concat(...Object.values(SYN.data.linksByQuery)).sort((a, b) => a.score - b.score);
    for (const l of allCross) {
        // skip chromosomes dropped by the min-length filter
        if (!(l.q_chrom in queryOffsets) || !(l.s_chrom in subjectOffsets)) { continue; }
        const aq1 = SYN.bpToAngleJS(l.q_chrom, l.q_start, SYN.data.querySizes, queryOffsets);
        const aq2 = SYN.bpToAngleJS(l.q_chrom, l.q_end, SYN.data.querySizes, queryOffsets);
        const as1 = SYN.bpToAngleJS(l.s_chrom, l.s_start, SYN.data.subjectSizes, subjectOffsets);
        const as2 = SYN.bpToAngleJS(l.s_chrom, l.s_end, SYN.data.subjectSizes, subjectOffsets);
        const poly = SYN.ribbonPolygonJS(aq1, aq2, ...SYN.orientEnds(as1, as2, l), SYN.data.linkR);
        ribbons.push({
            xs: poly.xs, ys: poly.ys, palette_index: SYN.data.subjectColorIndex[l.s_chrom],
            alpha: SYN.linkAlpha(l.score, SYN.data.maxScore), label: SYN.formatLinkLabel(l),
            score: l.score, is_homeolog: false, q_chrom: l.q_chrom, s_chrom: l.s_chrom, key: SYN.linkKey(l),
        });
    }
    const targetHomeologLinks = SYN.data.targetHomeologLinks.slice().sort((a, b) => a.score - b.score);
    for (const l of targetHomeologLinks) {
        if (!(l.q_chrom in queryOffsets) || !(l.s_chrom in queryOffsets)) { continue; }
        const aq1 = SYN.bpToAngleJS(l.q_chrom, l.q_start, SYN.data.querySizes, queryOffsets);
        const aq2 = SYN.bpToAngleJS(l.q_chrom, l.q_end, SYN.data.querySizes, queryOffsets);
        const as1 = SYN.bpToAngleJS(l.s_chrom, l.s_start, SYN.data.querySizes, queryOffsets);
        const as2 = SYN.bpToAngleJS(l.s_chrom, l.s_end, SYN.data.querySizes, queryOffsets);
        const poly = SYN.ribbonPolygonJS(aq1, aq2, ...SYN.orientEnds(as1, as2, l), SYN.data.linkR);
        ribbons.push({
            xs: poly.xs, ys: poly.ys, palette_index: null, genome: 'target',
            alpha: SYN.linkAlpha(l.score, SYN.data.maxScore),
            label: SYN.formatLinkLabel(l, 'Target', 'Target') + ' (homeolog)',
            score: l.score, is_homeolog: true, q_chrom: l.q_chrom, s_chrom: l.s_chrom, key: SYN.linkKey(l),
        });
    }
    const referenceHomeologLinks = SYN.data.referenceHomeologLinks.slice().sort((a, b) => a.score - b.score);
    for (const l of referenceHomeologLinks) {
        if (!(l.q_chrom in subjectOffsets) || !(l.s_chrom in subjectOffsets)) { continue; }
        const aq1 = SYN.bpToAngleJS(l.q_chrom, l.q_start, SYN.data.subjectSizes, subjectOffsets);
        const aq2 = SYN.bpToAngleJS(l.q_chrom, l.q_end, SYN.data.subjectSizes, subjectOffsets);
        const as1 = SYN.bpToAngleJS(l.s_chrom, l.s_start, SYN.data.subjectSizes, subjectOffsets);
        const as2 = SYN.bpToAngleJS(l.s_chrom, l.s_end, SYN.data.subjectSizes, subjectOffsets);
        const poly = SYN.ribbonPolygonJS(aq1, aq2, ...SYN.orientEnds(as1, as2, l), SYN.data.linkR);
        ribbons.push({
            xs: poly.xs, ys: poly.ys, palette_index: SYN.data.subjectColorIndex[l.q_chrom], genome: 'reference',
            alpha: SYN.linkAlpha(l.score, SYN.data.maxScore),
            label: SYN.formatLinkLabel(l, 'Reference', 'Reference') + ' (homeolog)',
            score: l.score, is_homeolog: true, q_chrom: l.q_chrom, s_chrom: l.s_chrom, key: SYN.linkKey(l),
        });
    }
    return ribbons;
};

// A chromosome's on-screen name: its plot name, or its original FASTA ID
// with "Original names" on (*OriginalNames lists only changed names).
SYN.displayName = function(genome, name) {
    if (!SYN.state.originalNames) { return name; }
    const names = genome === 'target' ? SYN.data.targetOriginalNames : SYN.data.referenceOriginalNames;
    return names[name] || name;
};

SYN.buildRingLayout = function(queryOrder, subjectOrder, k) {
    const groupGap = SYN.data.groupGap;
    const subjectOffsets = SYN.computeCircularOffsets(
        subjectOrder, SYN.data.subjectSizes, Math.PI - groupGap, groupGap);
    const queryOffsets = SYN.computeCircularOffsets(
        queryOrder, SYN.data.querySizes, Math.PI + groupGap, 2 * Math.PI - groupGap);

    const q = {xs: [], ys: [], fill_color: [], name: [], size_label: [], group: [], display_name: []};
    for (const name of queryOrder) {
        const [a0, a1] = queryOffsets[name];
        const poly = SYN.wedgePolygonJS(a0, a1, SYN.data.innerR, SYN.data.outerR);
        q.xs.push(poly.xs); q.ys.push(poly.ys);
        q.fill_color.push(SYN.data.targetGrey); q.name.push(name);
        q.size_label.push(`${(SYN.data.querySizes[name] / 1e6).toFixed(2)} Mb`);
        q.group.push('target'); q.display_name.push(SYN.displayName('target', name));
    }

    const s = {xs: [], ys: [], fill_color: [], name: [], size_label: [], group: [], palette_index: [],
               display_name: []};
    for (const name of subjectOrder) {
        const [a0, a1] = subjectOffsets[name];
        const poly = SYN.wedgePolygonJS(a0, a1, SYN.data.innerR, SYN.data.outerR);
        const pIdx = SYN.data.subjectColorIndex[name];
        s.xs.push(poly.xs); s.ys.push(poly.ys);
        s.fill_color.push(SYN.paletteColor(pIdx, k));
        s.name.push(name);
        s.size_label.push(`${(SYN.data.subjectSizes[name] / 1e6).toFixed(2)} Mb`);
        s.group.push('reference'); s.palette_index.push(pIdx);
        s.display_name.push(SYN.displayName('reference', name));
    }

    // wedges narrower than ringLabelMinAngle go unlabelled
    const labelR = SYN.data.outerR + 0.14;
    const label = {x: [], y: [], text: []};
    const addLabels = (order, offsets, genome) => {
        for (const name of order) {
            const [a0, a1] = offsets[name];
            if (Math.abs(a1 - a0) < SYN.data.ringLabelMinAngle) { continue; }
            const mid = (a0 + a1) / 2;
            label.x.push(labelR * Math.cos(mid)); label.y.push(labelR * Math.sin(mid));
            label.text.push(SYN.displayName(genome, name));
        }
    };
    addLabels(queryOrder, queryOffsets, 'target');
    addLabels(subjectOrder, subjectOffsets, 'reference');

    const ribbons = SYN.buildRibbonRecords(queryOffsets, subjectOffsets);

    // computed regardless of showGaps; SYN.applyGapVisibility filters
    const gaps = SYN.buildGapRecords(queryOffsets, subjectOffsets, SYN.GAP_ARC_SEGMENTS);

    return {q, s, label, ribbons, gaps, queryOffsets, subjectOffsets};
};

// A gap wedge is a fraction of a degree wide, so 2 arc segments look the same
// as wedgePolygonJS's default 48 and cut each wedge from 98 vertices to 8.
SYN.GAP_ARC_SEGMENTS = 2;

// Ring wedges for each genome's assembly gaps, full ideogram height,
// widened to SYN.data.minGapAngle.
SYN.buildGapRecords = function(queryOffsets, subjectOffsets, nSegments) {
    const gaps = [];
    for (const chrom in SYN.data.targetGapsByChrom) {
        if (!(chrom in queryOffsets)) { continue; }
        for (const g of SYN.data.targetGapsByChrom[chrom]) {
            const a0 = SYN.bpToAngleJS(chrom, g.start, SYN.data.querySizes, queryOffsets);
            const a1 = SYN.bpToAngleJS(chrom, g.end, SYN.data.querySizes, queryOffsets);
            const [lo, hi] = SYN.widenGapAngle(a0, a1);
            const poly = SYN.wedgePolygonJS(lo, hi, SYN.data.innerR, SYN.data.outerR, nSegments);
            gaps.push({xs: poly.xs, ys: poly.ys, label: SYN.formatGapLabel('target', chrom, g.start, g.end)});
        }
    }
    for (const chrom in SYN.data.referenceGapsByChrom) {
        if (!(chrom in subjectOffsets)) { continue; }
        for (const g of SYN.data.referenceGapsByChrom[chrom]) {
            const a0 = SYN.bpToAngleJS(chrom, g.start, SYN.data.subjectSizes, subjectOffsets);
            const a1 = SYN.bpToAngleJS(chrom, g.end, SYN.data.subjectSizes, subjectOffsets);
            const [lo, hi] = SYN.widenGapAngle(a0, a1);
            const poly = SYN.wedgePolygonJS(lo, hi, SYN.data.innerR, SYN.data.outerR, nSegments);
            gaps.push({xs: poly.xs, ys: poly.ys, label: SYN.formatGapLabel('reference', chrom, g.start, g.end)});
        }
    }
    return gaps;
};

// Applies a full ring reorder: rebuilds the wedges, labels and the ribbon
// and gap master lists, then re-applies the current filters and selection.
SYN.applyRingLayout = function(queryOrder, subjectOrder, sources) {
    const layout = SYN.buildRingLayout(queryOrder, subjectOrder, sources.k);
    sources.querySource.data = layout.q;
    sources.subjectSource.data = layout.s;
    sources.labelSource.data = layout.label;
    SYN.data.ribbons = layout.ribbons;
    SYN.data.gapRecords = layout.gaps;
    SYN.state.ringQueryOffsets = layout.queryOffsets;
    SYN.state.ringSubjectOffsets = layout.subjectOffsets;
    sources.querySource.change.emit();
    sources.subjectSource.change.emit();
    sources.labelSource.change.emit();
    SYN.applyOverviewRibbons(sources.minScore, sources.k, sources.showSelfLinks,
                              sources.showSynteny, sources.ribbonSource);
    SYN.applyGapVisibility(sources.gapSource);
};

// copies the current ring order's gap wedges into gap_src, or clears it
SYN.applyGapVisibility = function(gapSource) {
    const records = SYN.state.showGaps ? SYN.data.gapRecords : [];
    gapSource.data = {
        xs: records.map((r) => r.xs), ys: records.map((r) => r.ys), label: records.map((r) => r.label),
    };
    gapSource.change.emit();
};

// Dotplot block segments for the given axis offsets.
SYN.buildDotplotSegmentsForLayout = function(queryOffsets, subjectOffsets, minScore, k) {
    const xs = [], ys = [], lineColor = [], alpha = [], label = [], paletteIndex = [];
    const qChrom = [], sChrom = [], key = [];
    const maxScore = SYN.data.maxScore;
    for (const qName in SYN.data.linksByQuery) {
        for (const l of SYN.data.linksByQuery[qName]) {
            if (l.score < minScore) { continue; }
            // skip chromosomes dropped by the min-length filter
            if (!(l.q_chrom in queryOffsets) || !(l.s_chrom in subjectOffsets)) { continue; }
            const x0 = queryOffsets[l.q_chrom] + l.q_start;
            const x1 = queryOffsets[l.q_chrom] + l.q_end;
            const yLo = subjectOffsets[l.s_chrom] + l.s_start;
            const yHi = subjectOffsets[l.s_chrom] + l.s_end;
            // inverted blocks are drawn anti-diagonal
            const [y0, y1] = l.orientation !== '-' ? [yLo, yHi] : [yHi, yLo];
            xs.push([x0, x1]); ys.push([y0, y1]);
            const pIdx = SYN.data.subjectColorIndex[l.s_chrom];
            lineColor.push(SYN.linkColor(SYN.paletteColor(pIdx, k)));
            alpha.push(SYN.linkAlpha(l.score, maxScore));
            label.push(SYN.formatLinkLabel(l));
            paletteIndex.push(pIdx);
            qChrom.push(l.q_chrom); sChrom.push(l.s_chrom); key.push(SYN.linkKey(l));
        }
    }
    return {xs, ys, line_color: lineColor, alpha, label, palette_index: paletteIndex,
            q_chrom: qChrom, s_chrom: sChrom, key};
};

SYN.applyDotplotSegmentsForCurrentLayout = function(minScore, k, segmentSource) {
    segmentSource.data = SYN.buildDotplotSegmentsForLayout(SYN.state.dpQueryOffsets, SYN.state.dpSubjectOffsets, minScore, k);
    segmentSource.change.emit();
    SYN.reapplySelection();
};

SYN.computeOffsets = function(order, sizes) {
    const offsets = {};
    let cur = 0;
    for (const name of order) { offsets[name] = cur; cur += sizes[name]; }
    return offsets;
};

// The non-similarity order: size or FASTA order per SYN.state.orderBySize,
// first name leftmost/bottommost on both axes.
SYN.dpNaturalOrder = function() {
    const bySize = SYN.state.orderBySize;
    return {
        queryOrder: (bySize ? SYN.data.queryNames : SYN.data.queryNamesNatural).slice(),
        subjectOrder: (bySize ? SYN.data.subjectNames : SYN.data.subjectNamesNatural).slice(),
    };
};

// Drops chromosomes shorter than SYN.state.minSeqSize from an order.
SYN.filterBySize = function(order, sizes) {
    return order.filter((name) => sizes[name] >= SYN.state.minSeqSize);
};

// name -> position in `order` (SYN.data.subjectColorIndex). `order` must be
// the full, unfiltered list: the zoom panel can draw a partner below the
// min-length filter. With `% k` wrapping, chromosomes sharing a color are k
// apart in this order, so it's rebuilt when size/file order changes -- but
// not for similarity order, which re-runs on every re-chain and shouldn't
// reshuffle colors (so there, same-colored neighbours are possible).
SYN.buildColorIndex = function(order) {
    const index = {};
    order.forEach((name, i) => { index[name] = i; });
    return index;
};

// similarity order or SYN.dpNaturalOrder's size/natural pair, length-filtered
SYN.dpOrderFor = function(useSimilarity) {
    const order = useSimilarity ? SYN.computeSimilarityOrder() : SYN.dpNaturalOrder();
    return {
        queryOrder: SYN.filterBySize(order.queryOrder, SYN.data.querySizes),
        subjectOrder: SYN.filterBySize(order.subjectOrder, SYN.data.subjectSizes),
    };
};

// The one order both panels show: the ring's halves run left to right like
// the dotplot's axes, so similarity order also untangles the ring.
SYN.activeOrder = function() {
    return SYN.dpOrderFor(SYN.state.orderBySimilarity);
};

// Lays out the ring and the dotplot again for SYN.activeOrder.
SYN.applyOrder = function() {
    const ui = SYN.ui;
    const order = SYN.activeOrder();
    SYN.applyRingLayout(order.queryOrder, order.subjectOrder, {
        querySource: ui.querySource, subjectSource: ui.subjectSource, labelSource: ui.labelSource,
        ribbonSource: ui.ribbonSource, gapSource: ui.gapSource,
        minScore: ui.minBlockSpinner.value, k: ui.colorSpinner.value,
        showSelfLinks: ui.selfLinksToggle.active, showSynteny: ui.showSyntenyToggle.active,
    });
    SYN.applyDotplotOrder(order.queryOrder, order.subjectOrder, ui.colorSpinner.value, ui.minBlockSpinner.value, {
        query: ui.dpQuerySource, subject: ui.dpSubjectSource, grid: ui.dpGridSource,
        queryLabel: ui.dpQueryLabelSource, subjectLabel: ui.dpSubjectLabelSource,
        segment: ui.dpSegmentSource, gap: ui.dpGapSource, fig: ui.dotplotFig,
    });
};

// Reorders both axes so shared synteny reads as a diagonal: each chromosome
// goes at the score-weighted mean bp position of its blocks along the other
// axis (so chr3L and chr3R on one chr3 are ordered along it, not tied). The
// axes alternate a few times so they converge together. Only drawn blocks
// and chromosomes count; one with no visible blocks sorts to the end.
SYN.computeSimilarityOrder = function() {
    const natural = SYN.dpNaturalOrder();
    let queryOrder = SYN.filterBySize(natural.queryOrder, SYN.data.querySizes);
    let subjectOrder = SYN.filterBySize(natural.subjectOrder, SYN.data.subjectSizes);
    const minScore = SYN.ui.minBlockSpinner.value;

    // links: this axis's chromosome -> its blocks; other: {chrom, start, end}
    // field names for the other axis's side of each block
    const reorder = (names, otherOrder, otherSizes, links, other) => {
        const offsets = SYN.computeOffsets(otherOrder, otherSizes);
        return names
            .map((name, i) => {
                let totalW = 0, weightedPos = 0;
                for (const l of (links[name] || [])) {
                    const x0 = offsets[l[other.chrom]];
                    if (l.score < minScore || x0 === undefined) { continue; }
                    totalW += l.score;
                    weightedPos += l.score * (x0 + (l[other.start] + l[other.end]) / 2);
                }
                return {name, i, avgPos: totalW > 0 ? weightedPos / totalW : Infinity};
            })
            // index tiebreak: deterministic order among chromosomes with no signal
            .sort((a, b) => (a.avgPos - b.avgPos) || (a.i - b.i))
            .map((r) => r.name);
    };

    const sSide = {chrom: 's_chrom', start: 's_start', end: 's_end'};
    const qSide = {chrom: 'q_chrom', start: 'q_start', end: 'q_end'};
    for (let iter = 0; iter < 3; iter++) {
        queryOrder = reorder(queryOrder, subjectOrder, SYN.data.subjectSizes, SYN.data.linksByQuery, sSide);
        subjectOrder = reorder(subjectOrder, queryOrder, SYN.data.querySizes, SYN.data.linksByReference, qSide);
    }
    return {queryOrder, subjectOrder};
};

// Every order-dependent dotplot source for (queryOrder, subjectOrder).
// Totals are recomputed each call: the min-length filter changes them.
SYN.buildDotplotLayout = function(queryOrder, subjectOrder, k) {
    // `|| 1` keeps the math finite when the filter empties a whole side
    const totalX = queryOrder.reduce((sum, name) => sum + SYN.data.querySizes[name], 0) || 1;
    const totalY = subjectOrder.reduce((sum, name) => sum + SYN.data.subjectSizes[name], 0) || 1;
    const rx = totalX * SYN.data.dpRulerFrac, ry = totalY * SYN.data.dpRulerFrac;
    const queryOffsets = SYN.computeOffsets(queryOrder, SYN.data.querySizes);
    const subjectOffsets = SYN.computeOffsets(subjectOrder, SYN.data.subjectSizes);

    const qXs = [], qYs = [], qName = [], qSize = [], qGroup = [], qDisplay = [];
    for (const name of queryOrder) {
        const size = SYN.data.querySizes[name], x0 = queryOffsets[name];
        qXs.push([x0, x0 + size, x0 + size, x0]); qYs.push([-ry, -ry, 0, 0]);
        qName.push(name); qSize.push((size / 1e6).toFixed(2) + ' Mb'); qGroup.push('target');
        qDisplay.push(SYN.displayName('target', name));
    }
    const sXs = [], sYs = [], sName = [], sSize = [], sGroup = [], sIdx = [], sDisplay = [];
    for (const name of subjectOrder) {
        const size = SYN.data.subjectSizes[name], y0 = subjectOffsets[name];
        sXs.push([-rx, -rx, 0, 0]); sYs.push([y0, y0 + size, y0 + size, y0]);
        sName.push(name); sSize.push((size / 1e6).toFixed(2) + ' Mb'); sGroup.push('reference');
        sIdx.push(SYN.data.subjectColorIndex[name]); sDisplay.push(SYN.displayName('reference', name));
    }

    const gridXs = [], gridYs = [];
    for (const name of queryOrder) {
        const x0 = queryOffsets[name];
        gridXs.push([x0, x0]); gridYs.push([-ry, totalY]);
    }
    // the min-length filter can empty a whole side
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
        qLabelX.push(x0 + size / 2); qLabelY.push(-ry * 1.7); qLabelText.push(SYN.displayName('target', name));
    }
    const sLabelX = [], sLabelY = [], sLabelText = [];
    for (const name of subjectOrder) {
        const size = SYN.data.subjectSizes[name], y0 = subjectOffsets[name];
        sLabelX.push(-rx * 1.7); sLabelY.push(y0 + size / 2); sLabelText.push(SYN.displayName('reference', name));
    }

    return {
        queryOffsets, subjectOffsets, totalX, totalY, rx, ry,
        queryRuler: {xs: qXs, ys: qYs, fill_color: queryOrder.map(() => SYN.data.targetGrey),
                     name: qName, size_label: qSize, group: qGroup, display_name: qDisplay},
        subjectRuler: {xs: sXs, ys: sYs, fill_color: sIdx.map((idx) => SYN.paletteColor(idx, k)),
                       name: sName, size_label: sSize, group: sGroup, palette_index: sIdx,
                       display_name: sDisplay},
        grid: {xs: gridXs, ys: gridYs},
        queryLabels: {x: qLabelX, y: qLabelY, text: qLabelText},
        subjectLabels: {x: sLabelX, y: sLabelY, text: sLabelText},
    };
};

// Applies a dotplot reorder and records the new order/offsets in SYN.state
// for later redraws. sources.fig, if given, gets its ranges refit to the new
// totals (the min-length filter can shrink them).
SYN.applyDotplotOrder = function(queryOrder, subjectOrder, k, minScore, sources) {
    const layout = SYN.buildDotplotLayout(queryOrder, subjectOrder, k);
    SYN.state.dpQueryOffsets = layout.queryOffsets;
    SYN.state.dpSubjectOffsets = layout.subjectOffsets;
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
        // room for the right-aligned reference labels left of the ruler: *3
        // fits a 5-character name like "chr13" (*2.6 clips it), and longer
        // names (e.g. original accessions) get ~0.28 rx per character
        const longest = layout.subjectLabels.text.reduce((m, t) => Math.max(m, t.length), 0);
        const xMargin = Math.max(3, 1.7 + 0.28 * longest);
        sources.fig.x_range.start = -layout.rx * xMargin; sources.fig.x_range.end = layout.totalX * 1.02;
        sources.fig.y_range.start = -layout.ry * 3; sources.fig.y_range.end = layout.totalY * 1.02;
    }
    SYN.applyDotplotSegmentsForCurrentLayout(minScore, k, sources.segment);
    if (sources.gap) { SYN.applyDotplotGapVisibility(sources.gap); }
};

// Binary search for the name in `order` whose [offset, offset+size) contains
// pos (dotplot chromosomes are laid flush, so any pos in [0, total) hits
// one), or null.
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

// The dotplot cell under (x, y), by binary search per axis (a hidden
// polygon per cell would be ~49k for axolotl), or null outside the grid --
// the ruler strips and margins, where a tap falls through to the rulers.
SYN.dpCellAt = function(x, y) {
    if (x < 0 || y < 0 || x >= SYN.state.dpTotalX || y >= SYN.state.dpTotalY) { return null; }
    const targetName = SYN.chromAtOffset(SYN.state.dpQueryOrder, SYN.state.dpQueryOffsets, SYN.data.querySizes, x);
    const subjectName = SYN.chromAtOffset(SYN.state.dpSubjectOrder, SYN.state.dpSubjectOffsets, SYN.data.subjectSizes, y);
    if (targetName === null || subjectName === null) { return null; }
    return {targetName, subjectName};
};

// Dotplot gap lines for the current order, read from SYN.state so
// show_gaps_toggle's callback needn't know the order.
SYN.buildDotplotGapLines = function() {
    const rx = SYN.state.dpRx, ry = SYN.state.dpRy;
    const totalX = SYN.state.dpTotalX, totalY = SYN.state.dpTotalY;
    const xs = [], ys = [], label = [];
    for (const chrom in SYN.data.targetGapsByChrom) {
        // dropped by the min-length filter
        if (!(chrom in SYN.state.dpQueryOffsets)) { continue; }
        for (const g of SYN.data.targetGapsByChrom[chrom]) {
            const mid = SYN.state.dpQueryOffsets[chrom] + (g.start + g.end) / 2;
            xs.push([mid, mid]); ys.push([-ry, totalY]);
            label.push(SYN.formatGapLabel('target', chrom, g.start, g.end));
        }
    }
    for (const chrom in SYN.data.referenceGapsByChrom) {
        if (!(chrom in SYN.state.dpSubjectOffsets)) { continue; }
        for (const g of SYN.data.referenceGapsByChrom[chrom]) {
            const mid = SYN.state.dpSubjectOffsets[chrom] + (g.start + g.end) / 2;
            xs.push([-rx, totalX]); ys.push([mid, mid]);
            label.push(SYN.formatGapLabel('reference', chrom, g.start, g.end));
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

// per genome, since target and reference chromosomes can share names
SYN.gapsForChrom = function(role, name) {
    return (role === 'target' ? SYN.data.targetGapsByChrom : SYN.data.referenceGapsByChrom)[name] || [];
};

// two lines, since a Label doesn't wrap
SYN.EMPTY_DETAIL_MESSAGE = 'Click a chromosome on the ring or the dotplot,\nor a square in the dotplot, to see it here.';

SYN.plural = function(n, noun) {
    return `${n} ${noun}${n === 1 ? '' : 's'}`;
};

// empty: true means nothing selected: SYN.applyDetail hides the axis, grid
// and row labels and shows SYN.EMPTY_DETAIL_MESSAGE.
SYN.emptyDetail = function(title) {
    return {
        bars: {xs: [], ys: [], fill_color: []},
        ribbons: {xs: [], ys: [], fill_color: [], alpha: [], label: [], key: []},
        labels: {x: [], y: [], text: [], color: []},
        gaps: {xs: [], ys: [], label: []},
        x_range: [-1, 1], y_range: [-0.5, 1.5],
        title: title || 'Detail',
        empty: true,
    };
};

// pivotSide is 'subject' (reference) or 'query' (target), whichever was
// clicked: that chromosome gets its own row, its partners are packed on the
// other. Reference is always the top row.
SYN.buildDetailData = function(pivotSide, pivotName, k, minScore) {
    if (!pivotName) { return SYN.emptyDetail(); }
    const isPivotSubject = pivotSide === 'subject';
    const links = ((isPivotSubject ? SYN.data.linksByReference[pivotName] : SYN.data.linksByQuery[pivotName]) || [])
        .filter((l) => l.score >= minScore);
    const pivotSize = isPivotSubject ? SYN.data.subjectSizes[pivotName] : SYN.data.querySizes[pivotName];
    const pivotColor = isPivotSubject ? SYN.paletteColor(SYN.data.subjectColorIndex[pivotName], k) : SYN.data.targetGrey;

    const {barH, topY, botY} = SYN.data.detailRows;
    const pivotRowY = isPivotSubject ? topY : botY;
    const otherRowY = isPivotSubject ? botY : topY;
    const pivotRole = isPivotSubject ? 'reference' : 'target';
    const otherRole = isPivotSubject ? 'target' : 'reference';

    // full-height gap ticks for the bar starting at x0 on row rowY
    const gapXs = [], gapYs = [], gapLabel = [];
    const addGapTicks = (role, name, x0, rowY, refSpan) => {
        for (const g of SYN.gapsForChrom(role, name)) {
            const [gs, ge] = SYN.widenGapBp(g.start, g.end, refSpan);
            gapXs.push([x0 + gs, x0 + ge, x0 + ge, x0 + gs]);
            gapYs.push([rowY, rowY, rowY + barH, rowY + barH]);
            gapLabel.push(SYN.formatGapLabel(role, name, g.start, g.end));
        }
    };

    if (links.length === 0) {
        const d = SYN.emptyDetail(`${SYN.displayName(pivotRole, pivotName)} · no blocks`);
        d.empty = false;
        d.bars = {xs: [[0, pivotSize, pivotSize, 0]],
                  ys: [[pivotRowY, pivotRowY, pivotRowY + barH, pivotRowY + barH]],
                  fill_color: [pivotColor]};
        d.labels = {x: [pivotSize / 2], y: [pivotRowY + barH / 2], text: [SYN.displayName(pivotRole, pivotName)],
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
        labelText.push(SYN.displayName(otherRole, name));
        labelColor.push(SYN.textColorFor(color));
        addGapTicks(otherRole, name, x0, otherRowY, xMax);
    }
    barXs.push([0, pivotSize, pivotSize, 0]);
    barYs.push([pivotRowY, pivotRowY, pivotRowY + barH, pivotRowY + barH]);
    barFill.push(pivotColor);
    labelX.push(pivotSize / 2);
    labelY.push(pivotRowY + barH / 2);
    labelText.push(SYN.displayName(pivotRole, pivotName));
    labelColor.push(SYN.textColorFor(pivotColor));
    addGapTicks(pivotRole, pivotName, 0, pivotRowY, xMax);

    const otherEdge = isPivotSubject ? botY + barH : topY;
    const pivotEdge = isPivotSubject ? topY : botY + barH;
    const maxScore = links.reduce((m, l) => Math.max(m, l.score), 1);
    const ribXs = [], ribYs = [], ribFill = [], ribAlpha = [], ribLabel = [], ribKey = [];
    const sorted = links.slice().sort((a, b) => a.score - b.score);
    for (const l of sorted) {
        const oName = otherName(l);
        const ox0 = otherOffsets[oName] + otherStart(l);
        const ox1 = otherOffsets[oName] + otherEnd(l);
        const px0 = pivotStart(l), px1 = pivotEnd(l);
        ribXs.push([ox0, ox1, ...SYN.orientEnds(px1, px0, l)]);
        ribYs.push([otherEdge, otherEdge, pivotEdge, pivotEdge]);
        const refColor = isPivotSubject ? pivotColor : SYN.paletteColor(SYN.data.subjectColorIndex[oName], k);
        ribFill.push(SYN.linkColor(refColor));
        ribAlpha.push(SYN.linkAlpha(l.score, maxScore));
        ribLabel.push(SYN.formatLinkLabel(l));
        ribKey.push(SYN.linkKey(l));
    }

    return {
        bars: {xs: barXs, ys: barYs, fill_color: barFill},
        ribbons: {xs: ribXs, ys: ribYs, fill_color: ribFill, alpha: ribAlpha, label: ribLabel, key: ribKey},
        labels: {x: labelX, y: labelY, text: labelText, color: labelColor},
        gaps: {xs: gapXs, ys: gapYs, label: gapLabel},
        x_range: [-xMax * 0.03, xMax * 1.05],
        y_range: [botY - 0.35, topY + barH + 0.35],
        // the row labels name the genomes, keeping the title short
        title: `${SYN.displayName(pivotRole, pivotName)} · ${SYN.plural(links.length, 'block')} with ${SYN.plural(otherNames.length, 'chromosome')}`,
        empty: false,
    };
};

// One dotplot cell (target x reference chromosome): always two bars. Zero
// blocks is an expected outcome here -- "no synteny between this pair" is
// worth showing.
SYN.buildPairDetailData = function(targetName, subjectName, k, minScore) {
    const allLinks = SYN.data.linksByQuery[targetName] || [];
    const links = allLinks.filter((l) => l.s_chrom === subjectName && l.score >= minScore);
    const targetSize = SYN.data.querySizes[targetName];
    const subjectSize = SYN.data.subjectSizes[subjectName];
    const subjectColor = SYN.paletteColor(SYN.data.subjectColorIndex[subjectName], k);

    const {barH, topY, botY} = SYN.data.detailRows;
    const xMax = Math.max(targetSize, subjectSize);
    const barXs = [[0, subjectSize, subjectSize, 0], [0, targetSize, targetSize, 0]];
    const barYs = [[topY, topY, topY + barH, topY + barH], [botY, botY, botY + barH, botY + barH]];
    const barFill = [subjectColor, SYN.data.targetGrey];
    const labelX = [subjectSize / 2, targetSize / 2];
    const labelY = [topY + barH / 2, botY + barH / 2];
    const subjectLabel = SYN.displayName('reference', subjectName);
    const targetLabel = SYN.displayName('target', targetName);
    const labelText = [subjectLabel, targetLabel];
    const labelColor = [SYN.textColorFor(subjectColor), SYN.textColorFor(SYN.data.targetGrey)];
    const baseX = [-xMax * 0.03, xMax * 1.05];
    const baseY = [botY - 0.35, topY + barH + 0.35];

    // both bars start at x=0
    const gapXs = [], gapYs = [], gapLabel = [];
    const addGapTicks = (role, name, rowY) => {
        for (const g of SYN.gapsForChrom(role, name)) {
            const [gs, ge] = SYN.widenGapBp(g.start, g.end, xMax);
            gapXs.push([gs, ge, ge, gs]);
            gapYs.push([rowY, rowY, rowY + barH, rowY + barH]);
            gapLabel.push(SYN.formatGapLabel(role, name, g.start, g.end));
        }
    };
    addGapTicks('reference', subjectName, topY);
    addGapTicks('target', targetName, botY);
    const gaps = {xs: gapXs, ys: gapYs, label: gapLabel};

    if (links.length === 0) {
        return {
            bars: {xs: barXs, ys: barYs, fill_color: barFill},
            ribbons: {xs: [], ys: [], fill_color: [], alpha: [], label: [], key: []},
            labels: {x: labelX, y: labelY, text: labelText, color: labelColor},
            gaps,
            x_range: baseX, y_range: baseY,
            title: `${subjectLabel} × ${targetLabel} · no blocks`,
            empty: false,
        };
    }

    const maxScore = links.reduce((m, l) => Math.max(m, l.score), 1);
    const ribXs = [], ribYs = [], ribFill = [], ribAlpha = [], ribLabel = [], ribKey = [];
    const sorted = links.slice().sort((a, b) => a.score - b.score);
    for (const l of sorted) {
        ribXs.push([l.q_start, l.q_end, ...SYN.orientEnds(l.s_end, l.s_start, l)]);
        ribYs.push([botY + barH, botY + barH, topY, topY]);
        ribFill.push(SYN.linkColor(subjectColor));
        ribAlpha.push(SYN.linkAlpha(l.score, maxScore));
        ribLabel.push(SYN.formatLinkLabel(l));
        ribKey.push(SYN.linkKey(l));
    }

    return {
        bars: {xs: barXs, ys: barYs, fill_color: barFill},
        ribbons: {xs: ribXs, ys: ribYs, fill_color: ribFill, alpha: ribAlpha, label: ribLabel, key: ribKey},
        labels: {x: labelX, y: labelY, text: labelText, color: labelColor},
        gaps,
        x_range: baseX, y_range: baseY,
        // reference first, matching the rows
        title: `${subjectLabel} × ${targetLabel} · ${SYN.plural(links.length, 'block')}`,
        empty: false,
    };
};

// topLabel/bottomLabel default to reference over target; a homeolog link's
// two sides are both in one genome, so its callers pass that genome twice.
SYN.formatLinkLabel = function(l, topLabel, bottomLabel) {
    topLabel = topLabel || 'Reference';
    bottomLabel = bottomLabel || 'Target';
    const sName = SYN.escapeHtml(SYN.displayName(topLabel.toLowerCase(), l.s_chrom));
    const qName = SYN.escapeHtml(SYN.displayName(bottomLabel.toLowerCase(), l.q_chrom));
    const stats = `${l.score} anchors (distinct loci)  ·  orientation=${l.orientation}`
        + `<br>avg identity: ${(l.mean_identity * 100).toFixed(1)}%`
        + `  ·  density: ${l.anchor_density.toFixed(1)} anchors/Mb`;
    const coords = `<table style="border-collapse:collapse">`
        + `<tr><td style="padding-right:6px">${topLabel}</td>`
        + `<td style="text-align:right">${sName}:${l.s_start.toLocaleString()}-${l.s_end.toLocaleString()}</td></tr>`
        + `<tr><td style="padding-right:6px">${bottomLabel}</td>`
        + `<td style="text-align:right">${qName}:${l.q_start.toLocaleString()}-${l.q_end.toLocaleString()}</td></tr>`
        + `</table>`;
    return `${coords}<br>${stats}`;
};

SYN.formatGapLabel = function(genome, chrom, start, end) {
    return `Assembly gap<br>${SYN.escapeHtml(SYN.displayName(genome, chrom))}:${start.toLocaleString()}-${end.toLocaleString()}`
        + ` (${(end - start).toLocaleString()} bp)`;
};

// gapSource is optional; when given, gets d.gaps only if showGaps is on
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
    const ui = SYN.ui;
    for (const m of ui.detailAxisModels) { m.visible = !d.empty; }
    for (const m of ui.detailRowLabels) { m.visible = !d.empty; }
    ui.detailMessage.text = d.empty ? (d.message || SYN.EMPTY_DETAIL_MESSAGE) : '';
};

// ---- the one selection shared by all three panels ----
//
// SYN.state.selection is one of:
//   null
//   {kind: 'chrom', side: 'query'|'subject', name}   a ring wedge or dotplot band
//   {kind: 'pair', target, subject}                  a dotplot grid square
//   {kind: 'block', key, target, subject}            a synteny ribbon/segment
//   {kind: 'selfblock', key, genome, a, b}           a self-link ribbon (ring only)
// SYN.select is its only writer and shows it everywhere at once: source
// selections (Bokeh dims the rest), the dotplot highlight and the detail
// panel. Every source rebuild calls SYN.reapplySelection, so a selection
// survives re-chains, filters and reorders while it still exists.

// an index no source has: selecting it dims every glyph (an empty selection
// would draw them all normally)
SYN.DIM_ALL = [1000000000];

// the selection, or what's left of it once the data changed under it: a
// chromosome filtered out (min sequence length) drops the selection; a block
// that no longer exists after a re-chain or min-block change falls back to
// its chromosome pair; a self-link needs Show self-links on and to still exist
SYN.validSelection = function(sel) {
    if (!sel) { return null; }
    const ui = SYN.ui;
    const onRing = (side, name) => (side === 'query' ? ui.querySource : ui.subjectSource).data.name.includes(name);
    const minScore = ui.minBlockSpinner.value;
    if (sel.kind === 'chrom') { return onRing(sel.side, sel.name) ? sel : null; }
    if (!(sel.kind === 'selfblock' || (onRing('query', sel.target) && onRing('subject', sel.subject)))) { return null; }
    if (sel.kind === 'block') {
        const exists = (SYN.data.crossLinksFlat || []).some((l) => l.score >= minScore && SYN.linkKey(l) === sel.key);
        return exists ? sel : {kind: 'pair', target: sel.target, subject: sel.subject};
    }
    if (sel.kind === 'selfblock') {
        const side = sel.genome === 'target' ? 'query' : 'subject';
        const links = sel.genome === 'target' ? SYN.data.targetHomeologLinks : SYN.data.referenceHomeologLinks;
        const exists = ui.selfLinksToggle.active && onRing(side, sel.a) && onRing(side, sel.b)
            && links.some((l) => l.score >= minScore && SYN.linkKey(l) === sel.key);
        return exists ? sel : null;
    }
    return sel;
};

// what the detail panel shows for a selection
SYN.buildSelectionDetail = function(sel, k, minScore) {
    if (!sel) { return SYN.emptyDetail(); }
    if (sel.kind === 'chrom') { return SYN.buildDetailData(sel.side, sel.name, k, minScore); }
    if (sel.kind === 'selfblock') {
        const d = SYN.emptyDetail();
        const label = sel.genome === 'target' ? SYN.state.targetLabel : SYN.state.referenceLabel;
        d.message = `${SYN.displayName(sel.genome, sel.a)} × ${SYN.displayName(sel.genome, sel.b)}: a self-link within ${label}.\nSelf-links are drawn on the ring only.`;
        return d;
    }
    return SYN.buildPairDetailData(sel.target, sel.subject, k, minScore);
};

// [source, indices] for every selectable source
SYN.selectionIndices = function(sel) {
    const ui = SYN.ui;
    const byName = (src, names) => {
        const idx = [];
        src.data.name.forEach((n, i) => { if (names.has(n)) { idx.push(i); } });
        return idx.length ? idx : SYN.DIM_ALL;
    };
    const where = (src, pred) => {
        const idx = [];
        for (let i = 0; i < src.data.xs.length; i++) { if (pred(src.data, i)) { idx.push(i); } }
        return idx.length ? idx : SYN.DIM_ALL;
    };
    const none = (src) => [src, []];
    const all = [ui.querySource, ui.subjectSource, ui.ribbonSource, ui.selfRibbonSource,
                 ui.dpQuerySource, ui.dpSubjectSource, ui.dpSegmentSource, ui.detailRibbonSource];
    if (!sel) { return all.map(none); }
    const minScore = ui.minBlockSpinner.value;
    if (sel.kind === 'chrom') {
        const isSubject = sel.side === 'subject';
        const own = new Set([sel.name]);
        // the other genome's chromosomes this one has visible blocks with
        const links = ((isSubject ? SYN.data.linksByReference : SYN.data.linksByQuery)[sel.name] || [])
            .filter((l) => l.score >= minScore);
        const partners = new Set(links.map((l) => (isSubject ? l.q_chrom : l.s_chrom)));
        const col = isSubject ? 's_chrom' : 'q_chrom';
        const genome = isSubject ? 'reference' : 'target';
        return [
            [ui.querySource, byName(ui.querySource, isSubject ? partners : own)],
            [ui.subjectSource, byName(ui.subjectSource, isSubject ? own : partners)],
            [ui.ribbonSource, where(ui.ribbonSource, (d, i) => d[col][i] === sel.name)],
            [ui.selfRibbonSource, where(ui.selfRibbonSource, (d, i) => d.genome[i] === genome
                && (d.q_chrom[i] === sel.name || d.s_chrom[i] === sel.name))],
            [ui.dpQuerySource, byName(ui.dpQuerySource, isSubject ? partners : own)],
            [ui.dpSubjectSource, byName(ui.dpSubjectSource, isSubject ? own : partners)],
            [ui.dpSegmentSource, where(ui.dpSegmentSource, (d, i) => d[col][i] === sel.name)],
            none(ui.detailRibbonSource),
        ];
    }
    if (sel.kind === 'selfblock') {
        const isTarget = sel.genome === 'target';
        const pair = new Set([sel.a, sel.b]);
        return [
            [ui.querySource, isTarget ? byName(ui.querySource, pair) : SYN.DIM_ALL],
            [ui.subjectSource, isTarget ? SYN.DIM_ALL : byName(ui.subjectSource, pair)],
            [ui.ribbonSource, SYN.DIM_ALL],
            [ui.selfRibbonSource, where(ui.selfRibbonSource, (d, i) => d.key[i] === sel.key)],
            [ui.dpQuerySource, isTarget ? byName(ui.dpQuerySource, pair) : SYN.DIM_ALL],
            [ui.dpSubjectSource, isTarget ? SYN.DIM_ALL : byName(ui.dpSubjectSource, pair)],
            [ui.dpSegmentSource, SYN.DIM_ALL],
            none(ui.detailRibbonSource),
        ];
    }
    // pair or block
    const isBlock = sel.kind === 'block';
    const inPair = (d, i) => d.q_chrom[i] === sel.target && d.s_chrom[i] === sel.subject;
    const match = isBlock ? (d, i) => d.key[i] === sel.key : inPair;
    return [
        [ui.querySource, byName(ui.querySource, new Set([sel.target]))],
        [ui.subjectSource, byName(ui.subjectSource, new Set([sel.subject]))],
        [ui.ribbonSource, where(ui.ribbonSource, match)],
        [ui.selfRibbonSource, SYN.DIM_ALL],
        [ui.dpQuerySource, byName(ui.dpQuerySource, new Set([sel.target]))],
        [ui.dpSubjectSource, byName(ui.dpSubjectSource, new Set([sel.subject]))],
        [ui.dpSegmentSource, where(ui.dpSegmentSource, match)],
        [ui.detailRibbonSource, isBlock ? where(ui.detailRibbonSource, (d, i) => d.key[i] === sel.key) : []],
    ];
};

// the dotplot's highlight rectangle: a chromosome's whole row or column, or
// one grid square for a pair/block; none otherwise (or when the chromosome
// isn't on the dotplot's current axes)
SYN.selectionHighlight = function(sel) {
    const empty = {left: [], right: [], bottom: [], top: []};
    if (!sel || sel.kind === 'selfblock') { return empty; }
    const qo = SYN.state.dpQueryOffsets, so = SYN.state.dpSubjectOffsets;
    const span = (offsets, sizes, name) => (name in offsets ? [offsets[name], offsets[name] + sizes[name]] : null);
    let x = [0, SYN.state.dpTotalX], y = [0, SYN.state.dpTotalY];
    if (sel.kind === 'chrom') {
        if (sel.side === 'query') { x = span(qo, SYN.data.querySizes, sel.name); }
        else { y = span(so, SYN.data.subjectSizes, sel.name); }
    } else {
        x = span(qo, SYN.data.querySizes, sel.target);
        y = span(so, SYN.data.subjectSizes, sel.subject);
    }
    if (!x || !y) { return empty; }
    return {left: [x[0]], right: [x[1]], bottom: [y[0]], top: [y[1]]};
};

SYN.select = function(sel) {
    const ui = SYN.ui;
    sel = SYN.validSelection(sel);
    SYN.state.selection = sel;
    const d = SYN.buildSelectionDetail(sel, ui.colorSpinner.value, ui.minBlockSpinner.value);
    SYN.applyDetail(d, ui.detailBarSource, ui.detailRibbonSource, ui.detailLabelSource, ui.detailFig,
                     ui.detailGapSource);
    for (const [src, indices] of SYN.selectionIndices(sel)) {
        src.selected.indices = indices;
    }
    ui.dpHighlightSource.data = SYN.selectionHighlight(sel);
    ui.dpHighlightSource.change.emit();
};

// after any source a selection is drawn on was rebuilt
SYN.reapplySelection = function() {
    SYN.select(SYN.state.selection);
};

// ---- clicks. Every TapTool only hit-tests (behavior='inspect'), so Bokeh's
// own selection logic never fights SYN.select. A click is two events: the
// TapTool callback per renderer hit (SYN.onTapHit records it), and the
// figure's Tap event (SYN.resolveTap), which acts one tick later, after all
// hit callbacks have run.

// the selection a click on `source` at its hit indices means (the last
// index: the topmost glyph, since each source draws in increasing order)
SYN.selectionForHit = function(source, indices) {
    const ui = SYN.ui, i = indices[indices.length - 1], d = source.data;
    if (source === ui.subjectSource || source === ui.dpSubjectSource) {
        return {kind: 'chrom', side: 'subject', name: d.name[i]};
    }
    if (source === ui.querySource || source === ui.dpQuerySource) {
        return {kind: 'chrom', side: 'query', name: d.name[i]};
    }
    if (source === ui.selfRibbonSource) {
        return {kind: 'selfblock', key: d.key[i], genome: d.genome[i], a: d.q_chrom[i], b: d.s_chrom[i]};
    }
    if (source === ui.detailRibbonSource) {
        const l = (SYN.data.crossLinksFlat || []).find((l) => SYN.linkKey(l) === d.key[i]);
        return l ? {kind: 'block', key: d.key[i], target: l.q_chrom, subject: l.s_chrom} : null;
    }
    // ring ribbons, dotplot segments
    return {kind: 'block', key: d.key[i], target: d.q_chrom[i], subject: d.s_chrom[i]};
};

// a TapTool callback: remember what this click hit. When one click hits
// several renderers the last call wins -- build_page() lists each TapTool's
// renderers so the one that should win comes last.
SYN.onTapHit = function(panel, source) {
    const indices = source.inspected.indices;
    if (indices.length) { SYN.state.pendingTap[panel] = SYN.selectionForHit(source, indices); }
};

// a figure's Tap event: select what the click hit, or else do what a click
// on nothing means there -- `fallback` returns that selection, or undefined
// to leave the selection alone
SYN.resolveTap = function(panel, fallback) {
    setTimeout(function() {
        const hit = SYN.state.pendingTap[panel];
        delete SYN.state.pendingTap[panel];
        if (hit) { SYN.select(hit); return; }
        const sel = fallback();
        if (sel !== undefined) { SYN.select(sel); }
    }, 0);
};

// A figure's <svg> sits inside nested shadow roots, which querySelector can't
// pierce, alongside toolbar icon <svg>s -- so collect them all and keep the
// largest (SYN.findMainSvg) rather than rely on Bokeh's internal class names.
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

// Bokeh.index only holds root views, so a figure's view is found by walking
// child_views (a Map keyed by model) down from each root.
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

// Exports a figure's live <svg> (output_backend="svg") as a standalone file.
SYN.exportFigureAsSVG = function(figModel, filename) {
    const view = SYN.findViewByModelId(figModel.id);
    const svgEl = view && view.el ? SYN.findMainSvg(view.el) : null;
    if (!svgEl) {
        console.error('SYN.exportFigureAsSVG: could not find a rendered <svg> for this figure.');
        return;
    }
    // a clone, so xmlns (missing on an in-page fragment) and <desc> don't
    // touch the original
    const clone = svgEl.cloneNode(true);
    if (!clone.getAttribute('xmlns')) {
        clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    }
    // the blocks TSV's settings line (bin/chain.js's OUTPUTS header), as the
    // file's description: vector editors keep it, the drawing ignores it
    const settings = SYN.currentSettings();
    const desc = document.createElementNS('http://www.w3.org/2000/svg', 'desc');
    desc.textContent = SYNCHAIN.settingsLine(settings.params, settings.extra).replace(/^# |\n$/g, '');
    clone.insertBefore(desc, clone.firstChild);
    const svgText = new XMLSerializer().serializeToString(clone);
    SYN.downloadBlob(new Blob([svgText], {type: 'image/svg+xml;charset=utf-8'}), filename);
};

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

// PNG/JPEG: rasterizes the figure's SVG on a canvas at `scale` x screen
// size. The white fill keeps transparent areas from going black in a JPEG
// (toBlob's quality arg is ignored for PNG).
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

// A self-describing export filename from the (sanitized, free-text) labels,
// ending in SYN.settingsTag -- the only record of the settings a PNG/JPEG keeps.
SYN.buildExportFilename = function(panel, format) {
    const sanitize = (s) => (s || '').replace(/[^a-zA-Z0-9_-]+/g, '_').replace(/^_+|_+$/g, '');
    const t = sanitize(SYN.state.targetLabel) || 'target';
    const r = sanitize(SYN.state.referenceLabel) || 'reference';
    const ext = format === 'JPEG' ? 'jpg' : format.toLowerCase();
    return `${t}_vs_${r}_${panel}_${SYN.settingsTag()}.${ext}`;
};

// Sets the editable genome labels and refreshes the titles/annotations and
// the 'group' column (hover tooltip only) that aren't rebuilt on redraw.
SYN.applyLabels = function(targetLabel, referenceLabel, ctx) {
    SYN.state.targetLabel = targetLabel;
    SYN.state.referenceLabel = referenceLabel;
    ctx.overview.title.text = `${referenceLabel} (reference) vs ${targetLabel} (target)`;
    ctx.dotplotFig.title.text = `${referenceLabel} vs ${targetLabel}: whole-genome dotplot`;
    // in-panel genome names, which exports keep
    const ui = SYN.ui;
    ui.ringReferenceLabel.text = `Reference · ${referenceLabel}`;
    ui.ringTargetLabel.text = `Target · ${targetLabel}`;
    ui.dotplotXAxis.axis_label = `Target · ${targetLabel}`;
    ui.dotplotYAxis.axis_label = `Reference · ${referenceLabel}`;
    ui.detailRowLabels[0].text = referenceLabel;
    ui.detailRowLabels[1].text = targetLabel;
    for (const src of ctx.queryLikeSources) {
        src.data.group = src.data.group.map(() => targetLabel);
        src.change.emit();
    }
    for (const src of ctx.subjectLikeSources) {
        src.data.group = src.data.group.map(() => referenceLabel);
        src.change.emit();
    }
};

// truncates to 1 decimal place, so 99.96% never shows as 100.0%
SYN.trunc1 = function(x) {
    return Math.floor(x * 10) / 10;
};

SYN.buildStatsHtml = function(targetLabel, referenceLabel) {
    const s = SYN.data.alignmentStats;
    if (!s) { return ''; }
    const esc = SYN.escapeHtml;
    const pct = (n) => s.proteome_total ? SYN.trunc1(n / s.proteome_total * 100).toFixed(1) : '0.0';
    // each input's file name or accession (*_source), species in italics
    // below; long names are ellipsized, full name in the tooltip
    const clip = 'overflow:hidden;text-overflow:ellipsis;white-space:nowrap';
    const fg2 = SYN.data.theme.fg2;
    const nameCell = (source, species) => {
        const src = source ? `<div style="${clip}" title="${esc(source)}">${esc(source)}</div>` : '';
        const sp = species
            ? `<div style="${clip};color:${fg2};font-size:12px" title="${esc(species)}"><i>${esc(species)}</i></div>` : '';
        return src + sp || '&mdash;';
    };
    const td = 'padding:2px 6px 2px 0;vertical-align:top';
    const num = `${td};text-align:right;white-space:nowrap`;
    const row = (role, source, species, aligned, identity) => `
        <tr><td style="${td};color:${fg2}">${role}</td>
            <td style="${td}">${nameCell(source, species)}</td>
            <td style="${num}">${aligned.toLocaleString()} (${pct(aligned)}%)</td>
            <td style="${num};padding-right:0">${SYN.trunc1(identity * 100).toFixed(1)}%</td></tr>`;
    // width = STATS_PANEL_WIDTH, the zoom panel's drawn frame
    return `
        <div style="font-size:13px;color:${SYN.data.theme.fg};border:1px solid ${SYN.data.theme.line};
                    border-radius:12px;padding:10px 14px;display:block;box-sizing:border-box;
                    width:${SYN.data.statsPanelWidth}px;margin:6px 0 4px 0">
            <div style="font-size:11px;font-weight:700;letter-spacing:.06em;text-transform:uppercase;
                        color:${fg2}">Alignment summary</div>
            <table style="border-collapse:collapse;margin-top:4px;width:100%;table-layout:fixed">
                <colgroup><col style="width:66px"><col><col style="width:100px"><col style="width:74px"></colgroup>
                <tr><td style="${td};color:${fg2}">Proteome</td>
                    <td style="${td}">${nameCell(s.proteome_source || s.proteome_origin, s.proteome_species)}</td>
                    <td colspan="2" style="${num};padding-right:0">${s.proteome_total.toLocaleString()} proteins</td></tr>
                <tr style="color:${fg2}"><td></td><td></td>
                    <td style="${num}">aligned</td><td style="${num};padding-right:0">avg identity</td></tr>
                ${row('Target', s.query_source || targetLabel, s.query_species, s.query_aligned, s.query_mean_identity)}
                ${row('Reference', s.subject_source || referenceLabel, s.subject_species,
                      s.subject_aligned, s.subject_mean_identity)}
            </table>
        </div>`;
};

SYN.applyStats = function(targetLabel, referenceLabel, statsDiv) {
    statsDiv.text = SYN.buildStatsHtml(targetLabel, referenceLabel);
};

// ---- client-side chaining (bin/chain.js): decode the embedded hit tables
// and chain them (in a Web Worker when available) at load and on every
// min-identity/max-gap/hit-rank change ----

// the chainer always runs at this minBlock (bin/chain.js's CHAINING RULES,
// rule 6: extraction doesn't depend on it); the min block size control only
// filters the reply, never re-chains
SYN.CHAIN_REQUEST_MIN_BLOCK = 3;

SYN.HIT_RANK_MAP = {'best only': 1, '≤ 2': 2, '≤ 3': 3, all: SYNCHAIN.DEFAULTS.maxHitRank};

SYN.now = function() {
    return (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
};

// the chain-affecting controls (not min block size, see above)
SYN.currentChainParams = function() {
    return {
        minPositive: SYN.ui.minIdentitySpinner.value / 100,
        maxGap: SYN.ui.maxGapSpinner.value,
        maxHitRank: SYN.HIT_RANK_MAP[SYN.ui.hitRankSelect.value],
        minBlock: SYN.CHAIN_REQUEST_MIN_BLOCK,
    };
};

// Runs inside the Worker after chain.js and nothing else: no SYN.* there, so
// it mirrors SYN.buildChainers/SYN.runChainers (the main-thread fallback)
// with SYNCHAIN calls. ${...} is substituted when this string is built.
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
        self.postMessage({
            seq: msg.seq,
            withSelf: msg.withSelf,
            crossLinks: SYNCHAIN.blocksToLinks(cross),
            targetSelfLinks: msg.withSelf ? SYNCHAIN.blocksToLinks(chainers.targetSelf.run(p)) : [],
            referenceSelfLinks: msg.withSelf ? SYNCHAIN.blocksToLinks(chainers.referenceSelf.run(p)) : [],
            ms: performance.now() - t0,
        });
    }
};
`;

// Main-thread fallback only: cross (target vs reference) plus each genome's
// self-comparison (homeologs), the latter run only while Show self-links is
// on (SYN.dispatchChain).
SYN.buildChainers = function(target, reference) {
    return {
        cross: SYNCHAIN.createChainer(target, reference, {}),
        targetSelf: SYNCHAIN.createChainer(target, target, {selfMode: true}),
        referenceSelf: SYNCHAIN.createChainer(reference, reference, {selfMode: true}),
    };
};

SYN.runChainers = function(chainers, params, seq, withSelf) {
    const p = Object.assign({}, params, {minBlock: SYN.CHAIN_REQUEST_MIN_BLOCK});
    const t0 = SYN.now();
    const cross = chainers.cross.run(p);
    return {
        seq,
        withSelf,
        crossLinks: SYNCHAIN.blocksToLinks(cross),
        targetSelfLinks: withSelf ? SYNCHAIN.blocksToLinks(chainers.targetSelf.run(p)) : [],
        referenceSelfLinks: withSelf ? SYNCHAIN.blocksToLinks(chainers.referenceSelf.run(p)) : [],
        ms: SYN.now() - t0,
    };
};

// Decodes the payload and starts a Worker from the #synchain-src copy of
// chain.js plus SYN.WORKER_SHIM_SRC, falling back to chaining on the main
// thread if the browser refuses to create one (e.g. on file: pages).
SYN.startChainer = async function(payload) {
    const chroms = {
        target: {names: SYN.data.queryNamesNatural, sizes: SYN.data.queryNamesNatural.map((n) => SYN.data.querySizes[n])},
        reference: {names: SYN.data.subjectNamesNatural, sizes: SYN.data.subjectNamesNatural.map((n) => SYN.data.subjectSizes[n])},
    };
    const tables = await SYNCHAIN.decodePayload(payload, chroms);
    SYN.chain = {
        tables, seq: -1, appliedSeq: -1, busy: false,
        pendingSeq: null, pendingParams: null, worker: null, useWorker: true, mainThreadChainers: null,
        selfFresh: false,
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

// the empty detail panel's message doubles as a busy indicator, restored by
// SYN.applyChainResult's SYN.reapplySelection
SYN.showComputingStatus = function() {
    if (!SYN.state.selection) {
        SYN.ui.detailMessage.text = 'Computing synteny…';
    }
};

// At most one request in flight; while busy only the latest params are kept
// and sent when the reply arrives (SYN.onChainMessage). seq drops stale replies.
SYN.requestChain = function(params) {
    SYN.chain.seq++;
    const seq = SYN.chain.seq;
    if (SYN.chain.busy) {
        SYN.chain.pendingSeq = seq;
        SYN.chain.pendingParams = params;
        return;
    }
    SYN.dispatchChain(seq, params);
};

// Self-links only run while Show self-links is on, read at send time so a
// queued request picks up a switch flipped meanwhile.
SYN.dispatchChain = function(seq, params) {
    const withSelf = SYN.ui.selfLinksToggle.active;
    SYN.chain.busy = true;
    SYN.showComputingStatus();
    SYN.updateChainStatus();
    if (SYN.chain.useWorker) {
        SYN.chain.worker.postMessage({type: 'run', seq, params, withSelf});
    } else {
        SYN.onChainMessage(SYN.runChainers(SYN.chain.mainThreadChainers, params, seq, withSelf));
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
        SYN.dispatchChain(seq, params);
    } else {
        // nothing left in flight: clear a "finding self-links…" even when
        // this reply was stale and skipped above
        SYN.updateChainStatus();
    }
};

SYN.groupLinksBy = function(links, key) {
    const g = {};
    for (const l of links) {
        (g[l[key]] || (g[l[key]] = [])).push(l);
    }
    return g;
};

// The one place a chain result reaches SYN.data; re-renders every panel
// through the same apply* functions other controls use.
SYN.applyChainResult = function(msg) {
    const ui = SYN.ui;
    // kept flat, in chain.js's sorted order, for SYN.exportBlocksTsv
    SYN.data.crossLinksFlat = msg.crossLinks;
    SYN.data.linksByQuery = SYN.groupLinksBy(msg.crossLinks, 'q_chrom');
    SYN.data.linksByReference = SYN.groupLinksBy(msg.crossLinks, 's_chrom');
    SYN.data.targetHomeologLinks = msg.targetSelfLinks;
    SYN.data.referenceHomeologLinks = msg.referenceSelfLinks;
    // whether the self-links above match the current params (see
    // SYN.dispatchChain); Show self-links re-chains when they don't
    SYN.chain.selfFresh = msg.withSelf;
    // ribbon/segment opacity scales with the cross blocks only, so turning
    // self-links on never restyles the synteny already on screen
    let maxScore = 1;
    for (const l of msg.crossLinks) { if (l.score > maxScore) { maxScore = l.score; } }
    SYN.data.maxScore = maxScore;
    let maxAnyScore = maxScore;
    for (const l of msg.targetSelfLinks) { if (l.score > maxAnyScore) { maxAnyScore = l.score; } }
    for (const l of msg.referenceSelfLinks) { if (l.score > maxAnyScore) { maxAnyScore = l.score; } }

    if (SYN.state.orderBySimilarity) {
        // similarity order depends on the blocks, which just changed
        SYN.applyOrder();
    } else {
        // same order as before: only the ribbons and segments move
        SYN.data.ribbons = SYN.buildRibbonRecords(SYN.state.ringQueryOffsets, SYN.state.ringSubjectOffsets);
        SYN.applyOverviewRibbons(ui.minBlockSpinner.value, ui.colorSpinner.value, ui.selfLinksToggle.active,
                                  ui.showSyntenyToggle.active, ui.ribbonSource);
        SYN.applyDotplotSegmentsForCurrentLayout(ui.minBlockSpinner.value, ui.colorSpinner.value, ui.dpSegmentSource);
    }
    // high is the longest block. Bokeh clamps a value above high for good,
    // so the viewer's own value (minBlockWanted) is restored whenever high
    // allows, setting value/high in whichever order stays in bounds
    const spinner = ui.minBlockSpinner;
    const floor = spinner.low;
    const high = Math.max(floor, maxAnyScore);
    const value = Math.min(SYN.state.minBlockWanted != null ? SYN.state.minBlockWanted : spinner.value, high);
    // CustomJS callbacks run asynchronously, so min_block_callback tells
    // this value from a viewer's by the value itself, not by a flag
    if (value !== spinner.value) { SYN.state.minBlockFitted = value; }
    if (value <= spinner.high) { spinner.value = value; spinner.high = high; } else { spinner.high = high; spinner.value = value; }
    // nothing longer than the floor to filter: disable, and say why
    const stuck = maxAnyScore <= floor;
    spinner.disabled = stuck;
    const any = msg.crossLinks.length || msg.targetSelfLinks.length || msg.referenceSelfLinks.length;
    // after this update has finished: showing the note re-lays out the page,
    // and a repaint in the middle of it throws (BokehJS "Out of bounds")
    setTimeout(() => {
        ui.minBlockNote.text = !stuck ? '' : any
            ? `Fixed at ${floor}: no block is longer. Loosen Synteny detection for longer blocks.`
            : 'No blocks to filter at these settings.';
        ui.minBlockNote.visible = stuck;
    }, 0);
    // also replaces SYN.showComputingStatus's message
    SYN.reapplySelection();

    SYN.chain.lastMs = msg.ms;
    SYN.updateChainStatus();
};

// "N blocks" passing the min-block filter, plus the self-link count (or
// "none", or still chaining) while Show self-links is on; the last chain's
// time is in the hover title.
SYN.updateChainStatus = function() {
    const ui = SYN.ui;
    if (!SYN.data.crossLinksFlat || !SYN.chain) { return; }
    const minScore = ui.minBlockSpinner.value;
    const visible = SYN.data.crossLinksFlat.filter((l) => l.score >= minScore).length;
    let text = SYN.plural(visible, 'block');
    if (ui.selfLinksToggle.active) {
        if (SYN.chain.busy || !SYN.chain.selfFresh) {
            text += ' · finding self-links…';
        } else {
            const self = SYN.data.targetHomeologLinks.concat(SYN.data.referenceHomeologLinks)
                .filter((l) => l.score >= minScore).length;
            text += self ? ` · ${SYN.plural(self, 'self-link block')}` : ' · no self-links at these settings';
        }
    }
    const ms = `last computed in ${(SYN.chain.lastMs || 0).toFixed(0)} ms`;
    ui.chainStatusDiv.text = `<span title="${ms}">${text}</span>`;
};

// The cross-genome blocks on screen (min block size and min sequence length
// applied), in chain.js's sorted order. With the filters at their floors and
// "Original names" off, this is byte-identical to
// `node bin/chain_blocks.mjs` with the same settings.
SYN.exportBlocksTsv = function() {
    const minScore = SYN.ui.minBlockSpinner.value;
    const minLen = SYN.state.minSeqSize;
    const filtered = (SYN.data.crossLinksFlat || []).filter((l) => l.score >= minScore
        && SYN.data.querySizes[l.q_chrom] >= minLen && SYN.data.subjectSizes[l.s_chrom] >= minLen);
    // named as the page shows them
    const links = !SYN.state.originalNames ? filtered : filtered.map((l) => Object.assign({}, l, {
        q_chrom: SYN.displayName('target', l.q_chrom), s_chrom: SYN.displayName('reference', l.s_chrom),
    }));
    const settings = SYN.currentSettings();
    const tsv = SYNCHAIN.linksToTsv(links, settings.params, settings.extra);
    SYN.downloadBlob(new Blob([tsv], {type: 'text/tab-separated-values'}), SYN.buildBlocksExportFilename());
};

// The settings an export records: chain parameters (min block size at its
// control's value) and the page's filters once moved off their defaults.
SYN.currentSettings = function() {
    const params = Object.assign(SYN.currentChainParams(), {minBlock: SYN.ui.minBlockSpinner.value});
    const extra = [];
    // 10 kb: the min sequence length spinner's floor
    if (SYN.state.minSeqSize > 10000) { extra.push(['min_seq_size', SYN.state.minSeqSize]); }
    if (SYN.state.originalNames) { extra.push(['names', 'original']); }
    return {params, extra};
};

// the same settings, short enough for a filename: id70_gap25_min5, then
// _rank<N> below "all" and _minlen<X>Mb above the 10 kb floor
SYN.settingsTag = function() {
    const p = SYN.currentSettings().params;
    const rank = p.maxHitRank < SYNCHAIN.DEFAULTS.maxHitRank ? `_rank${p.maxHitRank}` : '';
    const minLen = SYN.state.minSeqSize > 10000 ? `_minlen${+(SYN.state.minSeqSize / 1e6).toFixed(3)}Mb` : '';
    return `id${Math.round(p.minPositive * 100)}_gap${p.maxGap}_min${p.minBlock}${rank}${minLen}`;
};

SYN.buildBlocksExportFilename = function() {
    const sanitize = (s) => (s || '').replace(/[^a-zA-Z0-9_-]+/g, '_').replace(/^_+|_+$/g, '');
    const t = sanitize(SYN.state.targetLabel) || 'target';
    const r = sanitize(SYN.state.referenceLabel) || 'reference';
    return `${t}_vs_${r}_blocks_${SYN.settingsTag()}.tsv`;
};

// The page's entry point, on DocumentReady: stashes `s` as SYN.ui, renders
// the ring/dotplot with zero links, seeds the min-identity/min-block controls
// (--min_identity/--min_block or SYNCHAIN.DEFAULTS), and requests the first chain.
SYN.init = function(s) {
    SYN.ui = s;
    SYN.state.targetLabel = SYN.data.targetLabelDefault;
    SYN.state.referenceLabel = SYN.data.referenceLabelDefault;
    if (s.statsDiv) { SYN.applyStats(SYN.state.targetLabel, SYN.state.referenceLabel, s.statsDiv); }

    // also shows the empty selection (SYN.reapplySelection, via the ring's
    // ribbon rebuild), which puts the detail panel in its empty state
    SYN.applyOrder();

    SYN.startChainer(SYN.hitsPayload).then(function() {
        const init = SYN.data.initialParams;
        const minPositive = init.minPositive != null ? init.minPositive : SYNCHAIN.DEFAULTS.minPositive;
        const minBlock = init.minBlock != null ? init.minBlock : SYNCHAIN.DEFAULTS.minBlock;
        // fires both spinners' callbacks -- harmless before the first chain
        s.minIdentitySpinner.value = Math.round(minPositive * 100);
        s.minBlockSpinner.value = minBlock;
        // no change event when it already shows minBlock: the viewer's value
        // to keep through re-chains starts here (see SYN.applyChainResult)
        SYN.state.minBlockWanted = minBlock;
        SYN.requestChain(SYN.currentChainParams());
    }).catch(function(e) {
        console.error('SYN.startChainer failed -- the page has no synteny to show:', e);
    });
};

"""


def build_page(ds, query_name, subject_name, query_subtitle=None, subject_subtitle=None,
               alignment_stats=None, hits_payload=None, chain_js_src='',
               min_identity=None, max_gap=25, min_block=None):
    # the reference falls back to "reference" rather than subject_name
    target_label_default, reference_label_default = default_labels(
        alignment_stats, query_subtitle, subject_subtitle, query_name, "reference")
    query_names = [n for n, _ in ds.query_chroms]
    subject_names = [n for n, _ in ds.subject_chroms]
    subject_index = {name: i for i, (name, _) in enumerate(ds.subject_chroms)}

    q_src, s_src, r_src, r_self_src, label_src = build_overview_sources()

    # ring gap wedges, filled by SYN.applyGapVisibility
    gap_src = ColumnDataSource(dict(xs=[], ys=[], label=[]))

    lim = OUTER_R + 0.22
    # no 'tap': it would add a second, unrestricted TapTool button next to
    # the hidden renderer-scoped one added below
    overview = figure(width=620, height=620, match_aspect=True,
                       x_range=Range1d(-lim, lim), y_range=Range1d(-lim, lim),
                       title=f"{reference_label_default} (reference) vs {target_label_default} (target)",
                       tools="pan,wheel_zoom,reset",
                       output_backend="svg")
    overview.axis.visible = False
    overview.grid.visible = False
    overview.toolbar.logo = None

    # self-links first, so the cross-genome synteny draws over them
    self_ribbon_renderer = overview.patches('xs', 'ys', source=r_self_src, fill_color='fill_color',
                                             fill_alpha='alpha', line_color='line_color', line_alpha=0.9,
                                             line_width=0.8, line_dash='dashed',
                                             nonselection_fill_alpha=RIBBON_DIM_ALPHA / 2,
                                             nonselection_line_alpha=RIBBON_DIM_ALPHA)
    ribbon_renderer = overview.patches('xs', 'ys', source=r_src, fill_color='fill_color',
                                        line_color=None, fill_alpha='alpha',
                                        nonselection_fill_alpha=RIBBON_DIM_ALPHA)
    wedge_dim = dict(nonselection_fill_alpha=WEDGE_DIM_ALPHA, nonselection_line_alpha=WEDGE_DIM_ALPHA)
    query_renderer = overview.patches('xs', 'ys', source=q_src, fill_color='fill_color',
                                       line_color='black', line_width=0.5, **wedge_dim)
    subject_renderer = overview.patches('xs', 'ys', source=s_src, fill_color='fill_color',
                                         line_color='black', line_width=0.5, **wedge_dim)
    overview.text('x', 'y', source=label_src, text_align='center', text_baseline='middle',
                  text_font_size='10px')
    # which half is which genome, as annotations so exports keep them; the
    # corners are well clear of the chromosome labels
    ring_group_label_style = dict(text_font_size='11px', text_font_style='bold', text_color='#333333')
    ring_reference_label = Label(x=-lim + 0.03, y=lim - 0.03, text_baseline='top',
                                 text=f"Reference · {reference_label_default}", **ring_group_label_style)
    ring_target_label = Label(x=-lim + 0.03, y=-lim + 0.03, text_baseline='bottom',
                              text=f"Target · {target_label_default}", **ring_group_label_style)
    overview.add_layout(ring_reference_label)
    overview.add_layout(ring_target_label)
    # drawn on top of the wedges; hover only, not in overview_tap, and it
    # doesn't block clicks to the wedge underneath
    gap_renderer = overview.patches('xs', 'ys', source=gap_src, fill_color='#000000',
                                     fill_alpha=0.55, line_color=None)

    # follow_mouse: a ribbon's centre can be outside a zoomed view
    overview.add_tools(HoverTool(renderers=[ribbon_renderer, self_ribbon_renderer], tooltips="@label{safe}",
                                  point_policy='follow_mouse'))
    overview.add_tools(HoverTool(renderers=[query_renderer, subject_renderer],
                                  tooltips=[("Chromosome", "@display_name"), ("Size", "@size_label"),
                                            ("Group", "@group")]))
    overview.add_tools(HoverTool(renderers=[gap_renderer], tooltips="@label{safe}",
                                  point_policy='follow_mouse'))
    # One hidden TapTool for every clickable renderer. A tool added via
    # add_tools() is never auto-activated, hence active_tap.
    # behavior='inspect': hit-test only (SYN.onTapHit); the last listed wins
    # a click that hits several -- cross ribbons draw over self-links, so
    # they win too
    overview_tap = TapTool(renderers=[self_ribbon_renderer, subject_renderer, query_renderer, ribbon_renderer],
                           behavior='inspect', visible=False,
                           callback=CustomJS(code="SYN.onTapHit('ring', cb_data.source);"))
    overview.add_tools(overview_tap)
    overview.toolbar.active_tap = overview_tap
    # double-click resets pan/zoom and clears the selection
    overview.js_on_event(DoubleTap, CustomJS(args=dict(fig=overview), code=f"""
        fig.x_range.start = {-lim}; fig.x_range.end = {lim};
        fig.y_range.start = {-lim}; fig.y_range.end = {lim};
        SYN.select(null);
    """))

    # Whole-genome dotplot: target along x, reference along y, a segment per
    # block. The ruler strips select chromosomes like the ring wedges do.
    (dp_q_src, dp_s_src, dp_seg_src, dp_grid_src, dp_q_label_src, dp_s_label_src) = build_dotplot_sources()

    # dotplot gap lines (SYN.applyDotplotGapVisibility)
    dp_gap_src = ColumnDataSource(dict(xs=[], ys=[], label=[]))

    # placeholder ranges, set by SYN.applyDotplotOrder before first paint
    dotplot_fig = figure(width=620, height=620,
                          x_range=Range1d(-1, 1), y_range=Range1d(-1, 1),
                          title=f"{reference_label_default} vs {target_label_default}: whole-genome dotplot",
                          tools="pan,wheel_zoom,reset",
                          output_backend="svg")
    # axes only for their titles: the rulers stand in for ticks
    dotplot_fig.xaxis.axis_label = f"Target · {target_label_default}"
    dotplot_fig.yaxis.axis_label = f"Reference · {reference_label_default}"
    dotplot_fig.axis.major_tick_line_color = None
    dotplot_fig.axis.minor_tick_line_color = None
    dotplot_fig.axis.major_label_text_font_size = '0px'
    dotplot_fig.axis.axis_line_color = None
    dotplot_fig.axis.axis_label_text_font_size = '11px'
    dotplot_fig.axis.axis_label_text_font_style = 'bold'
    dotplot_fig.axis.axis_label_text_color = '#333333'
    dotplot_fig.axis.axis_label_standoff = 0
    dotplot_fig.grid.visible = False
    dotplot_fig.toolbar.logo = None

    dotplot_fig.multi_line('xs', 'ys', source=dp_grid_src, line_color='black', line_width=1)
    # the selection's row/column strip or grid square (SYN.selectionHighlight),
    # under the segments so it frames them without covering them
    dp_highlight_src = ColumnDataSource(dict(left=[], right=[], bottom=[], top=[]))
    dotplot_fig.quad(left='left', right='right', bottom='bottom', top='top', source=dp_highlight_src,
                     fill_color='#000000', fill_alpha=0.12, line_color='#333333', line_width=1.5)
    dp_segment_renderer = dotplot_fig.multi_line('xs', 'ys', source=dp_seg_src, line_color='line_color',
                                                  line_alpha='alpha', line_width=2,
                                                  nonselection_line_alpha=RIBBON_DIM_ALPHA)
    # the same segments again, 8 px wide and invisible in every state: the
    # hover/click target, since a 2 px line is hard to hit with a mouse
    dp_segment_hit_renderer = dotplot_fig.multi_line('xs', 'ys', source=dp_seg_src, line_width=8,
                                                      line_alpha=0, selection_line_alpha=0,
                                                      nonselection_line_alpha=0)
    # subtle dotted lines, above the segments but below the ruler strips
    dp_gap_renderer = dotplot_fig.multi_line('xs', 'ys', source=dp_gap_src, line_color=GAP_LINE_COLOR,
                                              line_dash='dotted', line_width=0.75)
    dp_query_renderer = dotplot_fig.patches('xs', 'ys', source=dp_q_src, fill_color='fill_color',
                                             line_color='black', line_width=0.5, **wedge_dim)
    dp_subject_renderer = dotplot_fig.patches('xs', 'ys', source=dp_s_src, fill_color='fill_color',
                                               line_color='black', line_width=0.5, **wedge_dim)
    dotplot_fig.text('x', 'y', source=dp_q_label_src, text_align='center', text_baseline='top',
                      text_font_size='8px')
    dotplot_fig.text('x', 'y', source=dp_s_label_src, text_align='right', text_baseline='middle',
                      text_font_size='8px')

    dotplot_fig.add_tools(HoverTool(renderers=[dp_segment_hit_renderer], tooltips="@label{safe}",
                                     point_policy='follow_mouse'))
    dotplot_fig.add_tools(HoverTool(renderers=[dp_query_renderer, dp_subject_renderer],
                                     tooltips=[("Chromosome", "@display_name"), ("Size", "@size_label"),
                                               ("Group", "@group")]))
    # line_policy='interp': with 'nearest' the tooltip snaps to one of the
    # line's two endpoints
    dotplot_fig.add_tools(HoverTool(renderers=[dp_gap_renderer], tooltips="@label{safe}",
                                     point_policy='follow_mouse', line_policy='interp'))
    # as overview_tap; a grid-square tap is handled by the figure's Tap
    # event and SYN.dpCellAt instead
    dotplot_tap = TapTool(renderers=[dp_segment_hit_renderer, dp_query_renderer, dp_subject_renderer],
                          behavior='inspect', visible=False,
                          callback=CustomJS(code="SYN.onTapHit('dotplot', cb_data.source);"))
    dotplot_fig.add_tools(dotplot_tap)
    dotplot_fig.toolbar.active_tap = dotplot_tap
    # Double-click and Reset go to the current full extent: Bokeh's Reset
    # restores the Range1d(-1, 1) placeholder, and the min-length filter can
    # change the extent after load.
    dotplot_reset_callback = CustomJS(args=dict(fig=dotplot_fig), code="""
        fig.x_range.start = -SYN.state.dpRx * 3; fig.x_range.end = SYN.state.dpTotalX * 1.02;
        fig.y_range.start = -SYN.state.dpRy * 3; fig.y_range.end = SYN.state.dpTotalY * 1.02;
    """)
    dotplot_fig.js_on_event(DoubleTap, dotplot_reset_callback)
    # runs after Bokeh's own reset to the placeholder
    dotplot_fig.js_on_event(Reset, dotplot_reset_callback)

    detail_bar_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[]))
    detail_rib_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[], alpha=[], label=[], key=[]))
    detail_label_src = ColumnDataSource(dict(x=[], y=[], text=[], color=[]))
    detail_gap_src = ColumnDataSource(dict(xs=[], ys=[], label=[]))
    # Range1d placeholders (SYN.applyDetail sets them on every click), not
    # the default DataRange1d, whose async auto-fit would overwrite them
    detail_fig = figure(width=DETAIL_FIG_WIDTH, height=302, x_axis_label='position (Mb)',
                         min_border_left=DETAIL_FRAME_LEFT,
                         x_range=Range1d(-1, 1), y_range=Range1d(-1, 1),
                         title="Detail",
                         tools="pan,wheel_zoom,reset",
                         output_backend="svg")
    detail_fig.yaxis.visible = False
    # start as SYN.emptyDetail leaves them, so nothing flashes before SYN.init
    detail_fig.xaxis.visible = False
    detail_fig.grid.visible = False
    detail_fig.toolbar.logo = None
    # row labels: screen x so a pan never slides them away, data y at the
    # rows; annotations so exports keep them
    detail_row_label_style = dict(x=4, x_units='screen', text_font_size='10px', text_color='#555555',
                                  visible=False)
    detail_reference_row_label = Label(y=DETAIL_TOP_Y + DETAIL_BAR_H + 0.04, text_baseline='bottom',
                                       text=reference_label_default, **detail_row_label_style)
    detail_target_row_label = Label(y=DETAIL_BOT_Y - 0.04, text_baseline='top',
                                    text=target_label_default, **detail_row_label_style)
    detail_fig.add_layout(detail_reference_row_label)
    detail_fig.add_layout(detail_target_row_label)
    # the empty panel's message, centred in SYN.emptyDetail's range
    detail_message = Label(x=0, y=0.5, text='', text_align='center', text_baseline='middle',
                           text_font_size='12px', text_color='#888888', text_line_height=1.4)
    detail_fig.add_layout(detail_message)
    # coordinates are bp; ticks show Mb
    detail_fig.xaxis.formatter = CustomJSTickFormatter(code="""
        const mb = tick / 1e6;
        return (Math.round(mb * 10) / 10).toString();
    """)
    detail_ribbon_renderer = detail_fig.patches('xs', 'ys', source=detail_rib_src,
                                                 fill_color='fill_color', line_color=None,
                                                 fill_alpha='alpha', nonselection_fill_alpha=RIBBON_DIM_ALPHA)
    detail_bar_renderer = detail_fig.patches('xs', 'ys', source=detail_bar_src,
                                              fill_color='fill_color', line_color='black', line_width=0.5)
    # drawn over the bars; hover only, as on the ring
    detail_gap_renderer = detail_fig.patches('xs', 'ys', source=detail_gap_src, fill_color='#000000',
                                              fill_alpha=0.55, line_color=None)
    detail_fig.text('x', 'y', source=detail_label_src, text_align='center', text_baseline='middle',
                     text_font_size='9px', text_color='color')
    # follow_mouse: a ribbon's centre can be outside a zoomed view
    detail_fig.add_tools(HoverTool(renderers=[detail_ribbon_renderer], tooltips="@label{safe}",
                                    point_policy='follow_mouse'))
    detail_fig.add_tools(HoverTool(renderers=[detail_gap_renderer], tooltips="@label{safe}",
                                    point_policy='follow_mouse'))
    # clicking a ribbon selects that block everywhere; as overview_tap
    detail_tap = TapTool(renderers=[detail_ribbon_renderer], behavior='inspect', visible=False,
                         callback=CustomJS(code="SYN.onTapHit('detail', cb_data.source);"))
    detail_fig.add_tools(detail_tap)
    detail_fig.toolbar.active_tap = detail_tap
    # Double-click and Reset go to the range SYN.applyDetail last set, as
    # for the dotplot
    detail_reset_callback = CustomJS(args=dict(fig=detail_fig), code="""
        const r = SYN.state.detailRange;
        if (!r) { return; }
        fig.x_range.start = r.x0; fig.x_range.end = r.x1;
        fig.y_range.start = r.y0; fig.y_range.end = r.y1;
    """)
    detail_fig.js_on_event(DoubleTap, detail_reset_callback)
    detail_fig.js_on_event(Reset, detail_reset_callback)

    target_label_input = TextInput(title="Target label", value=target_label_default,
                                    width=TOP_CONTROL_WIDTH)
    reference_label_input = TextInput(title="Reference label", value=reference_label_default,
                                       width=TOP_CONTROL_WIDTH)

    palette_select = Select(title="Color palette", value=DEFAULT_PALETTE_NAME,
                             options=list(PALETTES.keys()), width=TOP_CONTROL_WIDTH)
    color_spinner = Spinner(title="Colors", low=1, high=MAX_COLORS,
                             step=1, value=DEFAULT_COLORS, width=TOP_CONTROL_WIDTH)
    # Chain parameters (bin/chain.js's PARAMETERS header): min identity, max
    # gap and hit rank re-chain; min block size only filters. SYN.init sets
    # the real min identity/min block values; these placeholders equal
    # SYNCHAIN.DEFAULTS so nothing visibly changes.
    # the "?" icon next to a titled input's label, as a wrapped box (a
    # plain-text tooltip is one unwrapped line that runs off the page)
    def help_tip(text, position='bottom'):
        return Tooltip(content=HTML(f"<div style='width:260px;white-space:normal;line-height:1.35'>"
                                    f"{text}</div>"), position=position)

    min_identity_spinner = Spinner(title="Min identity (%)", low=30, high=100, step=1,
                                    value=round((0.7 if min_identity is None else min_identity) * 100),
                                    width=TOP_CONTROL_WIDTH,
                                    description=help_tip(
                                        "A protein's hit counts only if at least this share of its aligned "
                                        "residues are identical or similar (miniprot's Positive score). Starts "
                                        "at 70%. Lower it for distant species; raise it to cut noise from "
                                        "paralogs. The alignment summary's avg identity is a guide."))
    max_gap_spinner = Spinner(title="Max gap (anchors)", low=1, high=100, step=1,
                               value=max_gap, width=TOP_CONTROL_WIDTH,
                               description=help_tip(
                                   "The most anchors a block may skip between two consecutive matched "
                                   "anchors, on either genome. Counted in anchors, not base pairs, so it means "
                                   "the same in a compact genome and a huge one. Larger values will join "
                                   "fragmented blocks; smaller values will split blocks at small rearrangements."))
    HIT_RANK_OPTIONS = ["best only", "≤ 2", "≤ 3", "all"]
    hit_rank_select = Select(title="Hit rank", value="all", options=HIT_RANK_OPTIONS, width=TOP_CONTROL_WIDTH,
                             description=help_tip(
                                 "Which of each protein's alignments to use: only its best hit, or also "
                                 "its 2nd/3rd best, and so on. 'best only' gives the cleanest one-to-one "
                                 "synteny; keep lower-ranked hits to see duplicated copies, such as "
                                 "paralogs or a polyploid's subgenomes."))
    # low=3: a 2-anchor chain is barely evidence. high is a placeholder;
    # SYN.applyChainResult sets it to the largest block.
    min_block_spinner = Spinner(title="Min block size", low=3, high=1000,
                                 step=1, value=min_block or 15, width=TOP_CONTROL_WIDTH,
                                 description=help_tip(
                                     "The minimum number of anchors for a syntenic block to be drawn."))
    # shown when there's nothing to filter (SYN.applyChainResult)
    min_block_note = Div(text="", visible=False, width=TOP_CONTROL_WIDTH, margin=(0, 5, 5, 5),
                         styles={'font-size': '11.5px', 'line-height': '1.35', 'color': THEME['fg2']})
    # "N blocks" status (SYN.updateChainStatus), telling a slow re-chain
    # apart from "nothing matched"
    chain_status_div = Div(text="", width=TOP_CONTROL_WIDTH, align='end', margin=(0, 5, 12, 5),
                           stylesheets=[f"span{{display:inline-block;padding:3px 10px;border-radius:999px;"
                                        f"background:{THEME['accent_tint']};color:{THEME['accent_tint_fg']};"
                                        f"font-size:12.5px;font-weight:600;white-space:nowrap;}}"])
    # Client-side length filter on top of the pipeline's --min_seq_size, in
    # Mb. Never below MIN_SEQ_SIZE_FLOOR_MB; 0.1 Mb steps (0.01 -> 0.1 ->
    # 0.2), typed values snapped to 0.1 and capped at the largest chromosome
    # by min_seq_size_callback. Both bounds live there, not in low/high:
    # Bokeh refuses an arrow step past `low` rather than clamping, and
    # silently rejects a typed value above `high`. low=0 lets the step to 0
    # through for the callback to lift to the floor.
    max_seq_size = max(list(ds.query_sizes.values()) + list(ds.subject_sizes.values()), default=0)
    min_seq_size_spinner = Spinner(title="Min sequence length (Mb)", low=0,
                                    step=0.1, value=MIN_SEQ_SIZE_FLOOR_MB, width=TOP_CONTROL_WIDTH,
                                    description=help_tip(
                                        "Hides chromosomes and scaffolds shorter than this, in Mb, from "
                                        "every panel and from the blocks TSV download. Typed values "
                                        "round to the nearest 0.1 Mb. Sequences below the "
                                        "pipeline's <span style='white-space:nowrap'>--min_seq_size</span> "
                                        "(default 10 kb) were already left out of the page."),
                                    # at most 3 decimals, no trailing zeros: 0.01, 0.1, 0.2, 1
                                    format='0[.][000]')
    reset_btn = Button(label="✕ Clear selection", button_type="default", width=150,
                        height=TOOLBAR_CONTROL_HEIGHT)
    SAVE_BUTTON_WIDTH = 70
    # the ring's save button follows its switches, so it gets their divider
    # and is widened by it
    SAVE_DIVIDER_PAD = 6
    SAVE_DIVIDER_CSS = f':host{{border-left:1px solid var(--divider-color);padding-left:{SAVE_DIVIDER_PAD}px;}}'
    save_ring_btn = Button(label="⬇ save", button_type="default", width=SAVE_BUTTON_WIDTH + SAVE_DIVIDER_PAD + 1,
                            height=TOOLBAR_CONTROL_HEIGHT, stylesheets=[SAVE_DIVIDER_CSS])
    save_zoom_btn = Button(label="⬇ save", button_type="default", width=SAVE_BUTTON_WIDTH,
                            height=TOOLBAR_CONTROL_HEIGHT)
    save_dotplot_btn = Button(label="⬇ save", button_type="default", width=SAVE_BUTTON_WIDTH,
                               height=TOOLBAR_CONTROL_HEIGHT)
    # on the dotplot row, the panel showing every cross-genome block
    export_blocks_btn = Button(label="⬇ blocks TSV", button_type="default", width=110,
                                height=TOOLBAR_CONTROL_HEIGHT)
    EXPORT_FORMATS = ["SVG", "PNG", "JPEG"]
    DEFAULT_EXPORT_FORMAT = "PNG"
    ring_format_sel = Select(options=EXPORT_FORMATS, value=DEFAULT_EXPORT_FORMAT, width=80,
                              height=TOOLBAR_CONTROL_HEIGHT)
    zoom_format_sel = Select(options=EXPORT_FORMATS, value=DEFAULT_EXPORT_FORMAT, width=80,
                              height=TOOLBAR_CONTROL_HEIGHT)
    dotplot_format_sel = Select(options=EXPORT_FORMATS, value=DEFAULT_EXPORT_FORMAT, width=80,
                                 height=TOOLBAR_CONTROL_HEIGHT)
    # the ring's and dotplot's chromosome order; "Size" matches Dataset's
    order_select = Select(title="Chromosome order", value="Size", options=list(ORDER_MODES),
                          width=TOP_CONTROL_WIDTH,
                          description=help_tip(
                              "Size: largest first. File order: as in each FASTA file. Similarity: "
                              "chromosomes that share blocks sit next to each other, which turns the "
                              "dotplot's synteny into a diagonal and untangles the ring. Each chromosome "
                              "is placed at the average position of its blocks along the other genome, "
                              "weighted by their anchor counts. Only what is drawn counts: blocks below "
                              "Min block size and sequences below Min sequence length are ignored, and "
                              "chromosomes with no visible blocks go to the end."))
    # every switch: a fixed 28px track (Bokeh's stretches with the label) and
    # a tighter gap. SWITCH_WIDTH fits "Show self-links" on one line.
    SWITCH_WIDTH = 125
    SWITCH_CSS = ':host{gap:4px;align-items:center;}.bk-body{flex:0 0 28px;width:28px;}'
    # plus a leading divider, for a switch following another control
    SWITCH_DIVIDER_CSS = SWITCH_CSS + ':host{border-left:1px solid var(--divider-color);padding-left:6px;}'
    # off by default: self-links can clutter the ring (e.g. a polyploid).
    # Under the ring, the only panel that draws them.
    self_links_toggle = Switch(label="Show self-links", active=False, width=SWITCH_WIDTH,
                                stylesheets=[SWITCH_DIVIDER_CSS])
    # hides the cross-genome ribbons, independently of self_links_toggle
    show_synteny_toggle = Switch(label="Show synteny", active=True, width=SWITCH_WIDTH,
                                  stylesheets=[SWITCH_CSS])
    # off by default: a draft assembly can have thousands of gaps. Covers
    # all three panels. align='end' plus a taller bottom margin centres it
    # on the titled inputs' boxes.
    show_gaps_toggle = Switch(label="Show gaps", active=False, width=SWITCH_WIDTH,
                               stylesheets=[SWITCH_CSS], align='end', margin=(5, 5, 11, 5))
    # disabled when no name was changed in either genome
    original_names_toggle = Switch(label="Original names", active=False, width=SWITCH_WIDTH,
                                    disabled=not (ds.query_original_names or ds.subject_original_names),
                                    stylesheets=[SWITCH_CSS], align='end', margin=(5, 5, 11, 5))

    save_ring_btn.js_on_click(CustomJS(args=dict(fig=overview, fmt=ring_format_sel), code="""
        SYN.exportFigure(fig, 'ring', fmt.value);
    """))
    save_zoom_btn.js_on_click(CustomJS(args=dict(fig=detail_fig, fmt=zoom_format_sel), code="""
        SYN.exportFigure(fig, 'detail', fmt.value);
    """))
    save_dotplot_btn.js_on_click(CustomJS(args=dict(fig=dotplot_fig, fmt=dotplot_format_sel), code="""
        SYN.exportFigure(fig, 'dotplot', fmt.value);
    """))
    export_blocks_btn.js_on_click(CustomJS(code="SYN.exportBlocksTsv();"))

    # the color index follows the size/file choice only (SYN.buildColorIndex)
    order_callback = CustomJS(code="""
        const [bySize, bySimilarity] = SYN.data.orderModes[cb_obj.value];
        if (bySize !== SYN.state.orderBySize) {
            SYN.data.subjectColorIndex = SYN.buildColorIndex(bySize ? SYN.data.subjectNames : SYN.data.subjectNamesNatural);
        }
        SYN.state.orderBySize = bySize;
        SYN.state.orderBySimilarity = bySimilarity;
        SYN.applyOrder();
    """)
    order_select.js_on_change('value', order_callback)

    min_seq_size_callback = CustomJS(code=f"""
        // snap to 0.1 Mb within [floor, largest chromosome]: 0.11 -> 0.1,
        // 0.30000000000000004 -> 0.3. Setting the value re-runs this callback.
        const floor = {MIN_SEQ_SIZE_FLOOR_MB};
        const top = Math.max(floor, Math.floor({max_seq_size / 1e6} * 10) / 10);
        const snapped = Math.min(top, Math.max(floor, Math.round(cb_obj.value * 10) / 10));
        if (snapped !== cb_obj.value) {{
            cb_obj.value = snapped;
            return;
        }}
        SYN.state.minSeqSize = Math.round(cb_obj.value * 1e6);  // whole bp, drops the float drift
        SYN.applyOrder();
    """)
    min_seq_size_spinner.js_on_change('value', min_seq_size_callback)

    # ---- each figure's Tap event resolves a click (SYN.resolveTap): select
    # what was hit, or what a click on nothing means in that panel:
    # - ring: clear the selection
    ring_tap_callback = CustomJS(code="SYN.resolveTap('ring', () => null);")
    overview.js_on_event(Tap, ring_tap_callback)
    # - dotplot: a grid square selects that chromosome pair, even with no
    #   blocks; the margins do nothing
    dotplot_pair_tap_callback = CustomJS(code="""
        const cell = SYN.dpCellAt(cb_obj.x, cb_obj.y);
        SYN.resolveTap('dotplot', () => (cell ? {kind: 'pair', target: cell.targetName, subject: cell.subjectName}
                                              : undefined));
    """)
    dotplot_fig.js_on_event(Tap, dotplot_pair_tap_callback)
    # - detail: step a selected block back out to its pair (the rest of what
    #   the panel shows); otherwise nothing
    detail_tap_callback = CustomJS(code="""
        SYN.resolveTap('detail', () => {
            const sel = SYN.state.selection;
            return sel && sel.kind === 'block' ? {kind: 'pair', target: sel.target, subject: sel.subject} : undefined;
        });
    """)
    detail_fig.js_on_event(Tap, detail_tap_callback)

    color_callback = CustomJS(args=dict(
        subject_source=s_src, dp_subject_source=dp_s_src,
        overview_ribbon_source=r_src, dotplot_segment_source=dp_seg_src, min_block_spinner=min_block_spinner,
        self_links_toggle=self_links_toggle, show_synteny_toggle=show_synteny_toggle,
    ), code="""
        // the rebuilds re-apply the selection, recoloring the detail panel too
        const k = cb_obj.value;
        SYN.recolorSubjectWedges(k, subject_source);
        SYN.recolorSubjectWedges(k, dp_subject_source);
        SYN.applyOverviewRibbons(min_block_spinner.value, k, self_links_toggle.active,
                                  show_synteny_toggle.active, overview_ribbon_source);
        SYN.applyDotplotSegmentsForCurrentLayout(min_block_spinner.value, k, dotplot_segment_source);
    """)
    color_spinner.js_on_change('value', color_callback)

    min_block_callback = CustomJS(args=dict(
        overview_ribbon_source=r_src, dotplot_segment_source=dp_seg_src, color_spinner=color_spinner,
        self_links_toggle=self_links_toggle, show_synteny_toggle=show_synteny_toggle,
    ), code="""
        // the viewer's own value, kept through re-chains (SYN.applyChainResult)
        if (cb_obj.value === SYN.state.minBlockFitted) {
            SYN.state.minBlockFitted = null;
        } else {
            SYN.state.minBlockWanted = cb_obj.value;
        }
        const minScore = cb_obj.value;
        const k = color_spinner.value;
        SYN.applyOverviewRibbons(minScore, k, self_links_toggle.active,
                                  show_synteny_toggle.active, overview_ribbon_source);
        if (SYN.state.orderBySimilarity) {
            // similarity order depends on the min block size
            SYN.applyOrder();
        } else {
            SYN.applyDotplotSegmentsForCurrentLayout(minScore, k, dotplot_segment_source);
        }
        SYN.updateChainStatus();
    """)
    min_block_spinner.js_on_change('value', min_block_callback)

    # these re-chain; SYN.ui holds every widget they need
    chain_params_callback = CustomJS(code="SYN.requestChain(SYN.currentChainParams());")
    min_identity_spinner.js_on_change('value', chain_params_callback)
    max_gap_spinner.js_on_change('value', chain_params_callback)
    hit_rank_select.js_on_change('value', chain_params_callback)

    self_links_toggle_callback = CustomJS(args=dict(
        overview_ribbon_source=r_src, color_spinner=color_spinner, min_block_spinner=min_block_spinner,
        show_synteny_toggle=show_synteny_toggle,
    ), code="""
        // self-links are only chained while this is on: if stale, re-chain
        // (the cross chainer's cache makes its share free)
        if (cb_obj.active && !SYN.chain.selfFresh) {
            SYN.requestChain(SYN.currentChainParams());
            return;
        }
        SYN.applyOverviewRibbons(min_block_spinner.value, color_spinner.value, cb_obj.active,
                                  show_synteny_toggle.active, overview_ribbon_source);
        SYN.updateChainStatus();
    """)
    self_links_toggle.js_on_change('active', self_links_toggle_callback)

    show_synteny_toggle_callback = CustomJS(args=dict(
        overview_ribbon_source=r_src, color_spinner=color_spinner, min_block_spinner=min_block_spinner,
        self_links_toggle=self_links_toggle,
    ), code="""
        SYN.applyOverviewRibbons(min_block_spinner.value, color_spinner.value, self_links_toggle.active,
                                  cb_obj.active, overview_ribbon_source);
    """)
    show_synteny_toggle.js_on_change('active', show_synteny_toggle_callback)

    # the detail panel's gap ticks show/hide via SYN.reapplySelection
    show_gaps_toggle_callback = CustomJS(args=dict(gap_source=gap_src, dp_gap_source=dp_gap_src), code="""
        SYN.state.showGaps = cb_obj.active;
        SYN.applyGapVisibility(gap_source);
        SYN.applyDotplotGapVisibility(dp_gap_source);
        SYN.reapplySelection();
    """)
    show_gaps_toggle.js_on_change('active', show_gaps_toggle_callback)

    # names are baked into every source, so re-lay out everything
    original_names_callback = CustomJS(code="""
        SYN.state.originalNames = cb_obj.active;
        SYN.applyOrder();
        SYN.reapplySelection();
    """)
    original_names_toggle.js_on_change('active', original_names_callback)

    # swaps SYN.data.palette, which every recolor reads; composes with
    # color_spinner's count
    palette_callback = CustomJS(args=dict(
        subject_source=s_src, dp_subject_source=dp_s_src,
        overview_ribbon_source=r_src, dotplot_segment_source=dp_seg_src,
        color_spinner=color_spinner, min_block_spinner=min_block_spinner,
        self_links_toggle=self_links_toggle, show_synteny_toggle=show_synteny_toggle,
    ), code="""
        SYN.data.palette = SYN.data.palettes[cb_obj.value];
        const k = color_spinner.value;
        SYN.recolorSubjectWedges(k, subject_source);
        SYN.recolorSubjectWedges(k, dp_subject_source);
        SYN.applyOverviewRibbons(min_block_spinner.value, k, self_links_toggle.active,
                                  show_synteny_toggle.active, overview_ribbon_source);
        SYN.applyDotplotSegmentsForCurrentLayout(min_block_spinner.value, k, dotplot_segment_source);
    """)
    palette_select.js_on_change('value', palette_callback)

    reset_callback = CustomJS(code="SYN.select(null);")
    reset_btn.js_on_click(reset_callback)

    # under each panel, only what acts on that panel alone: the ring's two
    # ribbon switches, then each panel's save button and format menu
    left_col = column(overview,
                       row(show_synteny_toggle, self_links_toggle, save_ring_btn, ring_format_sel))
    right_col = column(dotplot_fig,
                        row(save_dotplot_btn, dotplot_format_sel, export_blocks_btn))

    header_title_div = Div(text="<h1>quick_synteny</h1><p>Fast and annotation-free synteny "
                                 "visualisation using spliced protein alignments</p>",
                           margin=(0, 0, 14, 2),
                           stylesheets=[f"h1{{margin:0;font-size:24px;font-weight:650;letter-spacing:-0.01em;"
                                        f"color:{THEME['fg']};}}"
                                        f"p{{margin:4px 0 0;font-size:14px;color:{THEME['fg2']};}}"])

    # filled by SYN.applyStats (empty without --stats); width and left margin
    # line it up with the zoom panel's drawn frame
    stats_div = Div(text='', width=STATS_PANEL_WIDTH, margin=(0, 0, 0, DETAIL_FRAME_LEFT))

    # under the zoom panel, the shortest column
    mid_col = column(detail_fig, row(reset_btn, save_zoom_btn, zoom_format_sel), stats_div)

    label_callback = CustomJS(args=dict(
        target_label_input=target_label_input, reference_label_input=reference_label_input,
        overview=overview, dotplot_fig=dotplot_fig, stats_div=stats_div,
        q_src=q_src, dp_q_src=dp_q_src, s_src=s_src, dp_s_src=dp_s_src,
    ), code="""
        const targetLabel = target_label_input.value.trim() || 'target';
        const referenceLabel = reference_label_input.value.trim() || 'reference';
        SYN.applyLabels(targetLabel, referenceLabel, {
            overview, dotplotFig: dotplot_fig,
            queryLikeSources: [q_src, dp_q_src], subjectLikeSources: [s_src, dp_s_src],
        });
        SYN.reapplySelection();  // the detail panel's titles and messages name the genomes
        SYN.applyStats(targetLabel, referenceLabel, stats_div);
    """)
    target_label_input.js_on_change('value', label_callback)
    reference_label_input.js_on_change('value', label_callback)

    # a titled group of page-wide controls, with a leading divider except
    # first in its row; margin-left clears the "?" icon of the group to its
    # left, which can overflow a long title
    def control_group(title, *children, divider=True):
        css = (':host{border-left:1px solid var(--divider-color);padding-left:8px;margin-left:14px;}'
               if divider else '')
        header = Div(text=f"<span style='font-size:11px;font-weight:700;letter-spacing:.06em;"
                          f"text-transform:uppercase;color:{THEME['fg2']}'>{title}</span>",
                     margin=(4, 5, 0, 5))
        return column(header, row(*children), stylesheets=[css] if css else [])

    controls = column(
        row(control_group("Genomes", reference_label_input, target_label_input, divider=False),
            control_group("Synteny detection (recomputes)", min_identity_spinner, max_gap_spinner,
                          hit_rank_select, chain_status_div)),
        row(control_group("Filters", column(min_block_spinner, min_block_note), min_seq_size_spinner, divider=False),
            control_group("Display", palette_select, color_spinner, order_select, show_gaps_toggle,
                          original_names_toggle)),
        stylesheets=[CARD_CSS])
    for col in (left_col, mid_col, right_col):
        col.stylesheets = [*col.stylesheets, CARD_CSS]
    layout = column(
        header_title_div,
        controls,
        # row()'s default spacing is 0
        row(left_col, mid_col, right_col, spacing=10),
        spacing=12,
    )
    # the app's look on every element (see THEME_CSS): first, so a
    # widget's own stylesheets still override it
    theme = InlineStyleSheet(css=THEME_CSS)
    for m in layout.references():
        if isinstance(m, UIElement):
            m.stylesheets = [theme, *m.stylesheets]
    for fig in (overview, detail_fig, dotplot_fig):
        fig.title.text_font = UI_FONT
        fig.title.text_color = THEME['fg']
        fig.outline_line_color = THEME['line']

    # An explicit Document so SYN.init can run on DocumentReady, which fires
    # once Bokeh has built every model and before a viewer can touch a
    # control. Its args become SYN.ui, the registry later code reads through.
    doc = Document()
    doc.add_root(layout)
    doc.js_on_event(DocumentReady, CustomJS(args=dict(
        q_src=q_src, s_src=s_src, r_src=r_src, r_self_src=r_self_src, label_src=label_src, gap_src=gap_src,
        dp_q_src=dp_q_src, dp_s_src=dp_s_src, dp_grid_src=dp_grid_src,
        dp_q_label_src=dp_q_label_src, dp_s_label_src=dp_s_label_src,
        dp_seg_src=dp_seg_src, dp_gap_src=dp_gap_src, dotplot_fig=dotplot_fig,
        color_spinner=color_spinner, min_block_spinner=min_block_spinner, min_block_note=min_block_note,
        min_identity_spinner=min_identity_spinner, max_gap_spinner=max_gap_spinner,
        hit_rank_select=hit_rank_select, self_links_toggle=self_links_toggle,
        show_synteny_toggle=show_synteny_toggle, chain_status_div=chain_status_div,
        bar_source=detail_bar_src, ribbon_source=detail_rib_src, label_source=detail_label_src,
        detail_fig=detail_fig, detail_gap_source=detail_gap_src, stats_div=stats_div,
        ring_reference_label=ring_reference_label, ring_target_label=ring_target_label,
        dotplot_x_axis=dotplot_fig.xaxis[0], dotplot_y_axis=dotplot_fig.yaxis[0],
        # hidden while the detail panel is empty (see SYN.applyDetail)
        detail_axis_models=[detail_fig.xaxis[0], detail_fig.xgrid[0], detail_fig.ygrid[0]],
        detail_row_labels=[detail_reference_row_label, detail_target_row_label],
        detail_message=detail_message, dp_highlight_src=dp_highlight_src,
    ), code="""
        SYN.init({
            statsDiv: stats_div,
            querySource: q_src, subjectSource: s_src, ribbonSource: r_src, selfRibbonSource: r_self_src,
            labelSource: label_src, gapSource: gap_src,
            dpQuerySource: dp_q_src, dpSubjectSource: dp_s_src, dpGridSource: dp_grid_src,
            dpQueryLabelSource: dp_q_label_src, dpSubjectLabelSource: dp_s_label_src,
            dpSegmentSource: dp_seg_src, dpGapSource: dp_gap_src, dotplotFig: dotplot_fig,
            colorSpinner: color_spinner, minBlockSpinner: min_block_spinner, minBlockNote: min_block_note,
            minIdentitySpinner: min_identity_spinner, maxGapSpinner: max_gap_spinner,
            hitRankSelect: hit_rank_select, selfLinksToggle: self_links_toggle,
            showSyntenyToggle: show_synteny_toggle, chainStatusDiv: chain_status_div,
            detailBarSource: bar_source, detailRibbonSource: ribbon_source,
            detailLabelSource: label_source, detailFig: detail_fig, detailGapSource: detail_gap_source,
            ringReferenceLabel: ring_reference_label, ringTargetLabel: ring_target_label,
            dotplotXAxis: dotplot_x_axis, dotplotYAxis: dotplot_y_axis,
            detailAxisModels: detail_axis_models, detailRowLabels: detail_row_labels,
            detailMessage: detail_message, dpHighlightSource: dp_highlight_src,
        });
    """))

    # INLINE, not CDN, so the page renders offline
    page_html = file_html(doc, INLINE, title=f"{query_name} vs {subject_name} -- interactive synteny")
    # right after the opening <head>, which comes before any script: the
    # inlined BokehJS contains the literal text "</head>" too
    head, sep, tail = page_html.partition("<head>")
    page_html = head + sep + "\n" + PAGE_CSS + tail

    syn_data = {
        'theme': THEME,
        'palette': PALETTES[DEFAULT_PALETTE_NAME],
        'palettes': PALETTES,
        'targetGrey': TARGET_GREY,
        'targetSelfGrey': TARGET_SELF_GREY,
        # the label inputs' starting values
        'targetLabelDefault': target_label_default,
        'referenceLabelDefault': reference_label_default,
        'queryNames': query_names,
        'subjectNames': subject_names,
        # natural (FASTA) order, for "File order"; queryNames/subjectNames
        # are size order
        'queryNamesNatural': [n for n, _ in ds.query_chroms_natural],
        'subjectNamesNatural': [n for n, _ in ds.subject_chroms_natural],
        'querySizes': ds.query_sizes,
        'subjectSizes': ds.subject_sizes,
        'subjectColorIndex': subject_index,
        # empty until the first chain result (SYN.applyChainResult); SYN.init's
        # first render runs before it and reads these
        'linksByReference': {}, 'linksByQuery': {},
        'targetHomeologLinks': [], 'referenceHomeologLinks': [],
        'maxScore': 1,
        # one dict per genome, since target and reference chromosomes can
        # share names
        'targetGapsByChrom': group_gaps_by_chrom(ds.query_gaps),
        'referenceGapsByChrom': group_gaps_by_chrom(ds.subject_gaps),
        # renamed chromosome -> original sequence ID, for the "Original names"
        # switch -- only changed names of chromosomes in chrom.sizes, so a
        # fragmented assembly's filtered-out scaffolds add nothing
        'targetOriginalNames': ds.query_original_names,
        'referenceOriginalNames': ds.subject_original_names,
        # embedded so the JS geometry matches the figure extents set here
        'outerR': OUTER_R, 'innerR': INNER_R, 'linkR': LINK_R,
        'groupGap': GROUP_GAP, 'chromGap': CHROM_GAP, 'minGapAngle': MIN_GAP_ANGLE,
        'ringLabelMinAngle': RING_LABEL_MIN_ANGLE,
        'orderModes': ORDER_MODES,
        'detailFigWidth': DETAIL_FIG_WIDTH,
        'detailRows': {'barH': DETAIL_BAR_H, 'topY': DETAIL_TOP_Y, 'botY': DETAIL_BOT_Y},
        'statsPanelWidth': STATS_PANEL_WIDTH,
        'dpRulerFrac': DP_RULER_FRAC,
        # null without --stats: no stats panel
        'alignmentStats': alignment_stats,
        # --min_identity/--min_block, or null for SYNCHAIN.DEFAULTS; read once
        # by SYN.init. --max_gap always has a value.
        'initialParams': {'minPositive': min_identity, 'maxGap': max_gap, 'minBlock': min_block},
    }
    # bin/chain.js ahead of SHARED_JS. A second, inert text/plain copy lets
    # SYN.startChainer build the Worker's Blob from source, since a Worker
    # can't share a parsed <script> and fetch() fails on file:// pages.
    synchain_src_tag = f'<script type="text/plain" id="synchain-src">{chain_js_src}</script>\n'
    injected = (synchain_src_tag + "<script>\n" + chain_js_src + "\n" + SHARED_JS
                + "\nSYN.data = " + json.dumps(syn_data) + ";\nSYN.hitsPayload = "
                + json.dumps(hits_payload) + ";\n</script>\n</body>")
    # the last </body>: the inlined BokehJS contains the literal text too
    head, _, tail = page_html.rpartition("</body>")
    return head + injected + tail


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--query_name', required=True)
    parser.add_argument('--subject_name', required=True)
    parser.add_argument('--query_chrom_sizes', required=True)
    parser.add_argument('--subject_chrom_sizes', required=True)
    parser.add_argument('--target_hits', required=True,
                         help="bin/extract_hits.py output for the target genome "
                              "(format: that script's docstring) -- embedded as SYN.hitsPayload "
                              "and chained client-side; Python never chains anything")
    parser.add_argument('--reference_hits', required=True, help='ditto, for the reference genome')
    parser.add_argument('--min_identity', type=float, default=None,
                         help='initial Min identity (%%) control value, as a 0-1 fraction -- '
                              "default (unset): SYNCHAIN.DEFAULTS' 0.7")
    parser.add_argument('--max_gap', type=int, default=25,
                         help="initial Max gap (anchors) control value (bin/chain.js's PARAMETERS "
                              "header)")
    parser.add_argument('--min_block', type=int, default=None,
                         help='initial Min block size control value -- default (unset): '
                              "SYNCHAIN.DEFAULTS' 15")
    parser.add_argument('--stats', default=None,
                         help='optional compute_alignment_stats.py JSON -- shown as a small '
                              'always-visible summary panel on the page (proteome size, each '
                              "genome's aligned-protein count/mean identity) rather than baked "
                              'into any one exportable figure')
    parser.add_argument('--query_subtitle', default=None,
                         help="optional input file name/accession of the target genome -- "
                              "names it in the page's Target label when --stats has no species "
                              "for it (see default_labels)")
    parser.add_argument('--subject_subtitle', default=None,
                         help='ditto, under the subject/reference name')
    parser.add_argument('--target_gaps', default=None,
                         help="optional bin/rename_sequences.py --out_gaps output (chrom, start, "
                              "end TSV, no header) for the target genome -- drawn as gap markers "
                              'on the ring/zoom panel/dotplot when the page\'s "Show gaps" switch is on')
    parser.add_argument('--reference_gaps', default=None,
                         help='ditto, for the reference genome')
    parser.add_argument('--target_lookup', default=None,
                         help="optional bin/rename_sequences.py --out_lookup output for the target "
                              "genome -- each chromosome's original sequence ID, for the page's "
                              "\"Original names\" switch")
    parser.add_argument('--reference_lookup', default=None,
                         help='ditto, for the reference genome')
    parser.add_argument('--out_prefix', required=True)
    args = parser.parse_args()

    ds = Dataset(args.query_chrom_sizes, args.subject_chrom_sizes, args.target_gaps, args.reference_gaps,
                 args.target_lookup, args.reference_lookup)
    if not ds.query_chroms or not ds.subject_chroms:
        sys.exit("ERROR: no chromosomes to plot -- check the pipeline's --min_seq_size "
                  "isn't filtering out every sequence")

    alignment_stats = None
    if args.stats:
        with open(args.stats) as f:
            alignment_stats = json.load(f)

    hits_payload = build_hits_payload(args.target_hits, args.reference_hits,
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
    if args.target_gaps or args.reference_gaps:
        gap_bits.append(f"{len(ds.query_gaps)} target / {len(ds.subject_gaps)} reference assembly gap(s)")
    print(f"[plot_synteny_interactive] wrote {out_path} ({target_n} target / {reference_n} reference "
          f"hit(s) embedded, chained client-side"
          + ((f"; {'; '.join(gap_bits)}") if gap_bits else '') + ")",
          file=sys.stderr)


if __name__ == '__main__':
    main()
