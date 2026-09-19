#!/usr/bin/env python3
"""Render a two-genome synteny plot as a self-contained interactive HTML page
(Bokeh, no server -- all interactivity is CustomJS). This is the pipeline's
only plot output -- there is deliberately no static-image equivalent, so
that "explore, then export exactly the panel you want as an SVG" (see each
panel's save button, SYN.exportFigureAsSVG) can be the one workflow, rather
than needing separate static renders kept in sync with it.

Three-column layout: a full Circos-style ring on the left (every wedge/
ribbon individually hoverable), a linear zoom/detail panel in the middle, a
whole-genome dotplot on the right. Clicking a comparison-chromosome wedge, a
target-chromosome wedge, or a chromosome band on either dotplot axis zooms
the middle panel into that one chromosome's links -- whichever side you
didn't click gets packed (multiple chromosomes side by side if more than one
is involved), the clicked chromosome stays alone on its own row. Clicking a
single square in the dotplot grid instead zooms straight into that one
(target, comparison) chromosome pair -- always exactly two bars, no packing,
and this is also the only way to select a pair with zero links (see
SYN.buildPairDetailData), since the axis rulers can only ever select one
chromosome at a time. Comparison-genome chromosomes are always drawn on top
in every panel, regardless of which side triggered the zoom. A dropdown
picks which of a few curated color palettes (see PALETTES) comparison-genome
chromosomes cycle through, and a numeric spinner recolors every panel by
cycling through 1-MAX_COLORS discrete colors from whichever one is active
instead of one color per chromosome; a second spinner filters every panel
down to blocks with at least that many supporting anchors (floored at
MBA_SLIDER_MIN, the threshold the embedded links were actually generated at
-- see build_synteny.nf's *_FOR_SLIDER processes) -- see
build_overview_sources()'s ribbon_records and SYN.buildOverviewRibbons for
why this can be a pure client-side filter with no server and no second copy
of the chaining algorithm in JS. Spinners rather than sliders since the
filter is a plain score comparison that works for any integer, not just a
handful of steps -- a free-form numeric input doesn't imply a false ceiling
the way a slider's end-of-track does. Two text inputs let a viewer relabel
"target"/"comparison" to the actual genome/species names before exporting a
panel (see SYN.applyLabels) -- every exportable title/bar-label reads from
these rather than the pipeline's literal role tags. A stats panel (see
SYN.applyStats) shows the alignment summary --stats optionally provides,
since that no longer has anywhere to live inside an exported image the way
it did on the static plots this replaced.

The geometry/offset math for the zoom panel exists twice on purpose: once in
Python (build_overview_sources()/build_dotplot_sources()) and once as plain
JS (SHARED_JS), so the browser can rebuild an arbitrary (pivot side, pivot
chromosome, color count, min block anchors, chromosome order) combination
instead of only combinations precomputed ahead of time. Python only ships
the raw links (grouped by comparison-genome chromosome and by target
chromosome) plus chromosome sizes/colors, as embedded JSON.

This is not run through bin/'s usual container (the pygenomeviz biocontainer
this does run through has no Bokeh by default) -- see
modules/local/pygenomeviz_plot.nf for why this process pip-installs Bokeh
into it at runtime instead of using a dedicated image.
"""
import argparse
import html
import json
import math
import sys

from bokeh.embed import file_html
from bokeh.events import DoubleTap
from bokeh.layouts import column, row
from bokeh.models import (ColumnDataSource, HoverTool, TapTool, CustomJS, CustomJSTickFormatter,
                           Button, Div, Range1d, Select, Spinner, TextInput, Toggle)
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
# every option needs exactly MAX_COLORS entries -- SYN.paletteColor (see
# SHARED_JS) indexes with idx % k for any k up to MAX_COLORS, so a shorter
# list would read past its own end and recolor as undefined
PALETTES = {
    'Categorical (default)': PALETTE,
    'Colorblind-safe': PALETTE_COLORBLIND,
    'Pastel': PALETTE_PASTEL,
}
DEFAULT_PALETTE_NAME = 'Categorical (default)'
TARGET_GREY = '#999999'

OUTER_R = 1.0
RING_WIDTH = 0.045
INNER_R = OUTER_R - RING_WIDTH
LINK_R = INNER_R
GROUP_GAP = math.radians(4)
CHROM_GAP = math.radians(0.5)

MAX_COLORS = 10
DEFAULT_COLORS = MAX_COLORS

# fraction of each axis's total span given to the dotplot's chromosome-ruler
# strips (see build_dotplot_sources()) -- the clickable regions standing in
# for x/y axis ticks, since Bokeh has no native "clickable tick label"
DP_RULER_FRAC = 0.035

# The min-block-anchors spinner only works because --links is always
# generated at MBA_SLIDER_MIN (see modules/local/build_synteny.nf's
# *_FOR_SLIDER processes): find_blocks() in build_synteny_blocks.py discovers
# blocks largest-first per chromosome pair and stops once the next chain is
# smaller than its threshold, so a run at threshold T's output is exactly the
# score >= T subset of a run at any lower threshold's output. That means
# filtering the MBA_SLIDER_MIN links client-side by score reproduces exactly
# what a fresh run at any higher threshold would have produced, with no need
# to ship the actual chaining algorithm to the browser or precompute more
# than one file. There's no MAX/STEP alongside this: the spinner (see
# build_page()) accepts any integer up to this dataset's actual largest
# block, not a fixed set of steps.
MBA_SLIDER_MIN = 5


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


def read_links(path):
    """Header-driven: mean_identity/anchor_density (see build_synteny_blocks.py)
    come through automatically when present, and are None for older links.tsv
    files that predate those two columns -- no schema-version flag needed."""
    links = []
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        for line in f:
            row = dict(zip(header, line.rstrip('\n').split('\t')))
            links.append({
                'q_chrom': row['query_chrom'], 'q_start': int(row['query_start']),
                'q_end': int(row['query_end']), 's_chrom': row['subject_chrom'],
                's_start': int(row['subject_start']), 's_end': int(row['subject_end']),
                'score': int(row['score']), 'orientation': row['orientation'],
                'mean_identity': float(row['mean_identity']) if 'mean_identity' in row else None,
                'anchor_density': float(row['anchor_density']) if 'anchor_density' in row else None,
            })
    return links


def circular_offsets(chroms, angle_start, angle_end):
    """Each chromosome's (bp=0 angle, bp=size angle). angle_end may be less
    than angle_start -- this is why the subject (top) arc needs bp to follow
    DEcreasing angle for its coordinates to still read left-to-right on
    screen, while the query (bottom) arc already reads left-to-right with
    increasing angle."""
    n = len(chroms)
    direction = 1 if angle_end >= angle_start else -1
    total_gap = CHROM_GAP * max(n - 1, 0)
    available = max(abs(angle_end - angle_start) - total_gap, 0.01)
    total_bp = sum(s for _, s in chroms) or 1
    offsets = {}
    cur = angle_start
    for name, size in chroms:
        span = available * (size / total_bp)
        offsets[name] = (cur, cur + direction * span)
        cur += direction * (span + CHROM_GAP)
    return offsets


def linear_offsets(chroms, gap):
    """Each chromosome's cumulative bp=0 start position along a single axis,
    laid out end to end with a fixed gap between consecutive chromosomes --
    the dotplot panel's axes (see build_dotplot_figure()), as opposed to
    circular_offsets()'s angular layout for the ring. Returns (offsets dict,
    total span including trailing gaps)."""
    offsets = {}
    cur = 0
    for name, size in chroms:
        offsets[name] = cur
        cur += size + gap
    return offsets, max(cur - gap, 1)


def bp_to_angle(chrom, bp_pos, sizes, offsets):
    a0, a1 = offsets[chrom]
    size = sizes[chrom]
    frac = min(max(bp_pos / size, 0), 1) if size else 0
    return a0 + frac * (a1 - a0)


def point(angle, radius):
    return radius * math.cos(angle), radius * math.sin(angle)


def lighten(color, amount=0.35):
    """Blend a #rrggbb color toward white -- same math as SHARED_JS's
    SYN.lighten(), kept in sync deliberately (see module docstring)."""
    r, g, b = (int(color.lstrip('#')[i:i + 2], 16) for i in (0, 2, 4))
    blend = lambda c: round(c + (255 - c) * amount)
    return '#{:02x}{:02x}{:02x}'.format(blend(r), blend(g), blend(b))


def arc_points(a0, a1, radius, n=48):
    if a1 < a0:
        a0, a1 = a1, a0
    return [point(a0 + (a1 - a0) * i / n, radius) for i in range(n + 1)]


def wedge_polygon(a0, a1, r0=INNER_R, r1=OUTER_R, n=48):
    """Annular sector (a ring wedge) as a closed polygon: outer arc a0->a1,
    then inner arc back a1->a0, same shape matplotlib.patches.Wedge draws."""
    outer = arc_points(a0, a1, r1, n)
    inner = list(reversed(arc_points(a0, a1, r0, n)))
    poly = outer + inner
    return [p[0] for p in poly], [p[1] for p in poly]


def quad_bezier(p0, p1, p2, n=24):
    pts = []
    for i in range(n + 1):
        t = i / n
        x = (1 - t) ** 2 * p0[0] + 2 * (1 - t) * t * p1[0] + t ** 2 * p2[0]
        y = (1 - t) ** 2 * p0[1] + 2 * (1 - t) * t * p1[1] + t ** 2 * p2[1]
        pts.append((x, y))
    return pts


def ribbon_polygon(angle_q1, angle_q2, angle_s1, angle_s2, radius=LINK_R, n=24):
    """Straight chord across the query arc, quadratic Bezier (control at the
    ring's center) over to the subject arc, straight chord back, Bezier back
    to the start -- the bowtie ribbon crossing through the middle of the
    ring."""
    p_q1, p_q2 = point(angle_q1, radius), point(angle_q2, radius)
    p_s1, p_s2 = point(angle_s1, radius), point(angle_s2, radius)
    origin = (0.0, 0.0)
    poly = [p_q1, p_q2]
    poly += quad_bezier(p_q2, origin, p_s2, n)[1:]
    poly += [p_s1]
    poly += quad_bezier(p_s1, origin, p_q1, n)[1:]
    return [p[0] for p in poly], [p[1] for p in poly]


class Dataset:
    def __init__(self, query_chrom_sizes, subject_chrom_sizes, links_path,
                 target_homeolog_links_path=None, comparison_homeolog_links_path=None):
        self.query_chroms = read_chrom_sizes(query_chrom_sizes)
        self.subject_chroms = read_chrom_sizes(subject_chrom_sizes)
        self.query_sizes = dict(self.query_chroms)
        self.subject_sizes = dict(self.subject_chroms)

        # comparison (subject) occupies the upper half of the ring, target the
        # lower half. subject's angle args are given start > end so its bp
        # coordinates run left-to-right like the query arc's do, not
        # backwards -- see circular_offsets.
        self.subject_offsets = circular_offsets(self.subject_chroms, math.pi - GROUP_GAP, GROUP_GAP)
        self.query_offsets = circular_offsets(self.query_chroms, math.pi + GROUP_GAP, 2 * math.pi - GROUP_GAP)

        # linear (non-angular) layout for the dotplot panel -- target along
        # x, comparison along y, chromosomes flush end to end (no gap): the
        # boundary gridlines build_dotplot_sources() draws already mark each
        # edge, so a gap would only waste axis space, not add information
        self.dp_query_offsets, self.dp_query_span = linear_offsets(self.query_chroms, 0)
        # reversed relative to self.subject_chroms (and to the ring's subject
        # arrangement): linear_offsets gives the FIRST name in its input list
        # the lowest -- bottommost -- position, so without this the first
        # comparison chromosome would sit at the very bottom of the axis and
        # reading up the axis would run last-to-first -- reversing the input
        # here instead puts the first chromosome at the top, so the axis
        # reads top-to-bottom in the same natural (file) order the x-axis
        # already reads left-to-right. This is also this dataset's "natural"
        # dotplot order in the sense the "order by similarity" toggle uses
        # (see SHARED_JS's SYN.dpNaturalOrder, its client-side counterpart).
        self.dp_subject_offsets, self.dp_subject_span = linear_offsets(list(reversed(self.subject_chroms)), 0)

        self.query_colors = {name: TARGET_GREY for name, _ in self.query_chroms}
        self.subject_colors = {name: PALETTE[i % DEFAULT_COLORS]
                                for i, (name, _) in enumerate(self.subject_chroms)}

        self.links = [l for l in read_links(links_path)
                      if l['q_chrom'] in self.query_offsets and l['s_chrom'] in self.subject_offsets]
        self.target_homeolog_links = []
        if target_homeolog_links_path:
            self.target_homeolog_links = [l for l in read_links(target_homeolog_links_path)
                                           if l['q_chrom'] in self.query_offsets and l['s_chrom'] in self.query_offsets]
        self.comparison_homeolog_links = []
        if comparison_homeolog_links_path:
            self.comparison_homeolog_links = [l for l in read_links(comparison_homeolog_links_path)
                                               if l['q_chrom'] in self.subject_offsets and l['s_chrom'] in self.subject_offsets]

    def link_angles(self, link, q_offsets, q_sizes, s_offsets, s_sizes):
        return (
            bp_to_angle(link['q_chrom'], link['q_start'], q_sizes, q_offsets),
            bp_to_angle(link['q_chrom'], link['q_end'], q_sizes, q_offsets),
            bp_to_angle(link['s_chrom'], link['s_start'], s_sizes, s_offsets),
            bp_to_angle(link['s_chrom'], link['s_end'], s_sizes, s_offsets),
        )

    def max_score(self):
        return max((l['score'] for l in self.links + self.target_homeolog_links
                     + self.comparison_homeolog_links), default=1) or 1


def format_link_label(link, extra=''):
    """mean_identity/anchor_density are None for links.tsv files predating
    those columns -- omit that line rather than printing 'None'."""
    stats = f"{link['score']} protein alignment(s)  ·  orientation={link['orientation']}"
    if link.get('mean_identity') is not None and link.get('anchor_density') is not None:
        stats += (f"<br>avg identity: {link['mean_identity'] * 100:.1f}%"
                  f"  ·  density: {link['anchor_density']:.1f} anchors/Mb")
    return (f"{link['q_chrom']}:{link['q_start']:,}-{link['q_end']:,}"
            f" ↔ {link['s_chrom']}:{link['s_start']:,}-{link['s_end']:,}"
            f"<br>{stats}{extra}")


def format_stats_html(target_label, comparison_label, stats):
    """The interactive HTML's replacement for the (removed) static plots'
    stats inset -- see compute_alignment_stats.py. Mirrors SHARED_JS's
    SYN.buildStatsHtml exactly (same markup), since that's what re-renders
    this same panel client-side whenever a viewer edits the target/comparison
    labels -- kept in sync deliberately (see module docstring)."""
    if not stats:
        return ''
    total = stats['proteome_total']

    def row(label, aligned, identity):
        pct = (aligned / total * 100) if total else 0.0
        return (f"<tr><td>{html.escape(label)}</td><td>{aligned:,} ({pct:.1f}%)</td>"
                f"<td>{identity * 100:.1f}%</td></tr>")

    return (
        "<div style='font-size:13px;color:#333;border:1px solid #ddd;border-radius:6px;"
        "padding:8px 12px;display:inline-block;margin:4px 0 8px 0'>"
        f"<b>Alignment summary</b> -- {total:,} protein(s) in the input proteome"
        "<table style='border-collapse:collapse;margin-top:4px'>"
        "<tr style='color:#666'><th style='text-align:left;padding-right:12px'></th>"
        "<th style='text-align:right;padding-right:12px'>aligned</th>"
        "<th style='text-align:right'>avg identity</th></tr>"
        + row(target_label, stats['query_aligned'], stats['query_mean_identity'])
        + row(comparison_label, stats['subject_aligned'], stats['subject_mean_identity'])
        + "</table></div>"
    )


def build_overview_sources(ds, subject_index, query_index):
    q_xs, q_ys, q_fill, q_name, q_size, q_group = [], [], [], [], [], []
    for name, size in ds.query_chroms:
        a0, a1 = ds.query_offsets[name]
        xs, ys = wedge_polygon(a0, a1)
        q_xs.append(xs); q_ys.append(ys)
        q_fill.append(ds.query_colors[name]); q_name.append(name)
        q_size.append(f"{size / 1e6:.2f} Mb"); q_group.append('target')

    s_xs, s_ys, s_fill, s_name, s_size, s_group, s_idx = [], [], [], [], [], [], []
    for name, size in ds.subject_chroms:
        a0, a1 = ds.subject_offsets[name]
        xs, ys = wedge_polygon(a0, a1)
        s_xs.append(xs); s_ys.append(ys)
        s_fill.append(ds.subject_colors[name]); s_name.append(name)
        s_size.append(f"{size / 1e6:.2f} Mb"); s_group.append('comparison')
        s_idx.append(subject_index[name])

    # Every ribbon's full geometry/label is computed once here regardless of
    # score, and carried into the page as SYN.data.ribbons (see build_page())
    # so the min-block-anchors slider can filter by score and recolor
    # entirely client-side (SYN.buildOverviewRibbons) -- no data is dropped
    # at this stage.
    ribbon_records = []
    max_score = ds.max_score()
    for link in sorted(ds.links, key=lambda l: l['score']):
        aq1, aq2, as1, as2 = ds.link_angles(link, ds.query_offsets, ds.query_sizes,
                                             ds.subject_offsets, ds.subject_sizes)
        xs, ys = ribbon_polygon(aq1, aq2, as1, as2)
        # dotplot segment endpoints: same block, plotted in linear (target-x,
        # comparison-y) coordinates instead of angular ones. A '-' (inverted)
        # block's target and comparison spans run opposite directions, so the
        # segment is drawn corner-to-corner the other way (anti-diagonal) --
        # the same convention MCScanX/D-GENIES-style dotplots use to make
        # inversions visually distinct from collinear blocks.
        dp_x0 = ds.dp_query_offsets[link['q_chrom']] + link['q_start']
        dp_x1 = ds.dp_query_offsets[link['q_chrom']] + link['q_end']
        dp_y_lo = ds.dp_subject_offsets[link['s_chrom']] + link['s_start']
        dp_y_hi = ds.dp_subject_offsets[link['s_chrom']] + link['s_end']
        dp_y0, dp_y1 = (dp_y_lo, dp_y_hi) if link['orientation'] != '-' else (dp_y_hi, dp_y_lo)
        ribbon_records.append({
            'xs': xs, 'ys': ys,
            'palette_index': subject_index[link['s_chrom']],
            'alpha': 0.25 + 0.55 * (link['score'] / max_score),
            'label': format_link_label(link),
            'score': link['score'],
            'is_homeolog': False,
            'dp_x0': dp_x0, 'dp_y0': dp_y0, 'dp_x1': dp_x1, 'dp_y1': dp_y1,
        })
    # a genome's own self-comparison homeolog links (either genome, or both --
    # see Dataset.__init__) always land entirely within that genome's own
    # half of the ring, using its own offsets/sizes on both ends and its own
    # chromosome's palette index (query_index/subject_index) so the ribbon is
    # colored the same way that genome's wedges/ribbons already are. Neither
    # set has a natural place on a target-by-comparison dotplot -- excluded
    # there (the initial segment source below, and every later dotplot
    # redraw, which reads raw links from linksByQuery instead -- see
    # SHARED_JS's SYN.buildDotplotSegmentsForLayout -- rather than this
    # is_homeolog flag), so their dp_* fields are simply unused.
    for link in sorted(ds.target_homeolog_links, key=lambda l: l['score']):
        aq1, aq2, as1, as2 = ds.link_angles(link, ds.query_offsets, ds.query_sizes,
                                             ds.query_offsets, ds.query_sizes)
        xs, ys = ribbon_polygon(aq1, aq2, as1, as2)
        ribbon_records.append({
            'xs': xs, 'ys': ys,
            'palette_index': query_index[link['q_chrom']],
            'alpha': 0.25 + 0.55 * (link['score'] / max_score),
            'label': format_link_label(link, extra=' (homeolog)'),
            'score': link['score'],
            'is_homeolog': True,
            'dp_x0': 0, 'dp_y0': 0, 'dp_x1': 0, 'dp_y1': 0,
        })
    for link in sorted(ds.comparison_homeolog_links, key=lambda l: l['score']):
        aq1, aq2, as1, as2 = ds.link_angles(link, ds.subject_offsets, ds.subject_sizes,
                                             ds.subject_offsets, ds.subject_sizes)
        xs, ys = ribbon_polygon(aq1, aq2, as1, as2)
        ribbon_records.append({
            'xs': xs, 'ys': ys,
            'palette_index': subject_index[link['q_chrom']],
            'alpha': 0.25 + 0.55 * (link['score'] / max_score),
            'label': format_link_label(link, extra=' (homeolog)'),
            'score': link['score'],
            'is_homeolog': True,
            'dp_x0': 0, 'dp_y0': 0, 'dp_x1': 0, 'dp_y1': 0,
        })

    # query and subject chromosomes can share names (e.g. a polyploid comparison
    # renamed to chr1..chrN same as the target) -- iterate each side against
    # its own offsets dict rather than a single by-name lookup, which would
    # silently resolve every shared name to whichever dict happens to have it
    label_x, label_y, label_text = [], [], []
    for name, _size in ds.query_chroms:
        a0, a1 = ds.query_offsets[name]
        mid = (a0 + a1) / 2
        label_x.append((OUTER_R + 0.14) * math.cos(mid))
        label_y.append((OUTER_R + 0.14) * math.sin(mid))
        label_text.append(name)
    for name, _size in ds.subject_chroms:
        a0, a1 = ds.subject_offsets[name]
        mid = (a0 + a1) / 2
        label_x.append((OUTER_R + 0.14) * math.cos(mid))
        label_y.append((OUTER_R + 0.14) * math.sin(mid))
        label_text.append(name)

    # the overview starts showing every ribbon (slider defaults to
    # MBA_SLIDER_MIN, i.e. no filtering), so the initial ColumnDataSource is
    # just ribbon_records unfiltered, colored at the default color count
    r_src = ColumnDataSource(dict(
        xs=[r['xs'] for r in ribbon_records],
        ys=[r['ys'] for r in ribbon_records],
        fill_color=[lighten(PALETTE[r['palette_index'] % DEFAULT_COLORS]) for r in ribbon_records],
        alpha=[r['alpha'] for r in ribbon_records],
        label=[r['label'] for r in ribbon_records],
        palette_index=[r['palette_index'] for r in ribbon_records],
    ))

    return (
        ColumnDataSource(dict(xs=q_xs, ys=q_ys, fill_color=q_fill, name=q_name, size_label=q_size, group=q_group)),
        ColumnDataSource(dict(xs=s_xs, ys=s_ys, fill_color=s_fill, name=s_name, size_label=s_size,
                               group=s_group, palette_index=s_idx)),
        r_src,
        ColumnDataSource(dict(x=label_x, y=label_y, text=label_text)),
        ribbon_records,
    )


def build_dotplot_sources(ds, subject_index, ribbon_records):
    """Chromosome-ruler strips (clickable -- carry the same 'name' field the
    ring wedges use, so build_page() can wire them into the exact same
    zoom-pivot tap callback) plus boundary gridlines, the initial
    (unfiltered) block-segment source, and the invisible N x M grid-cell hit
    layer (one per (target chromosome, comparison chromosome) pair -- see the
    'cell_src' built below -- for the pairwise zoom, including pairs with no
    links at all), all in the dotplot's linear target-x / comparison-y
    coordinate space -- see Dataset.__init__'s dp_query_offsets/
    dp_subject_offsets and DP_RULER_FRAC. Homeolog links (either genome's own
    self-comparison) have no natural x-position on a target-by-comparison plot
    and are excluded here."""
    total_x, total_y = ds.dp_query_span, ds.dp_subject_span
    rx, ry = total_x * DP_RULER_FRAC, total_y * DP_RULER_FRAC

    # bottom ruler: one rect per target chromosome, just below y=0
    q_xs, q_ys, q_fill, q_name, q_size, q_group = [], [], [], [], [], []
    for name, size in ds.query_chroms:
        x0 = ds.dp_query_offsets[name]
        q_xs.append([x0, x0 + size, x0 + size, x0])
        q_ys.append([-ry, -ry, 0, 0])
        q_fill.append(ds.query_colors[name]); q_name.append(name)
        q_size.append(f"{size / 1e6:.2f} Mb"); q_group.append('target')

    # left ruler: one rect per comparison chromosome, just left of x=0
    s_xs, s_ys, s_fill, s_name, s_size, s_group, s_idx = [], [], [], [], [], [], []
    for name, size in ds.subject_chroms:
        y0 = ds.dp_subject_offsets[name]
        s_xs.append([-rx, -rx, 0, 0])
        s_ys.append([y0, y0 + size, y0 + size, y0])
        s_fill.append(ds.subject_colors[name]); s_name.append(name)
        s_size.append(f"{size / 1e6:.2f} Mb"); s_group.append('comparison')
        s_idx.append(subject_index[name])

    # faint boundary lines at every chromosome edge, spanning the full plot
    grid_xs, grid_ys = [], []
    for name, size in ds.query_chroms:
        x0 = ds.dp_query_offsets[name]
        grid_xs.append([x0, x0]); grid_ys.append([-ry, total_y])
    last_q_name, last_q_size = ds.query_chroms[-1]
    grid_xs.append([ds.dp_query_offsets[last_q_name] + last_q_size] * 2); grid_ys.append([-ry, total_y])
    for name, size in ds.subject_chroms:
        y0 = ds.dp_subject_offsets[name]
        grid_xs.append([-rx, total_x]); grid_ys.append([y0, y0])
    # the topmost edge belongs to ds.subject_chroms[0] here, not [-1] --
    # dp_subject_offsets was built from the reversed list (see
    # Dataset.__init__), so the FIRST name in the original list is the one
    # that ends up with the highest (topmost) position
    top_s_name, top_s_size = ds.subject_chroms[0]
    grid_xs.append([-rx, total_x]); grid_ys.append([ds.dp_subject_offsets[top_s_name] + top_s_size] * 2)

    q_label_x, q_label_y, q_label_text = [], [], []
    for name, size in ds.query_chroms:
        x0 = ds.dp_query_offsets[name]
        q_label_x.append(x0 + size / 2); q_label_y.append(-ry * 1.7); q_label_text.append(name)
    s_label_x, s_label_y, s_label_text = [], [], []
    for name, size in ds.subject_chroms:
        y0 = ds.dp_subject_offsets[name]
        s_label_x.append(-rx * 1.7); s_label_y.append(y0 + size / 2); s_label_text.append(name)

    # score >= MBA_SLIDER_MIN by construction (see Dataset), so the initial
    # segment source shows everything -- same convention as build_overview_
    # sources()'s r_src, kept in sync client-side by the same spinners
    cross = [r for r in ribbon_records if not r['is_homeolog']]
    seg_src = ColumnDataSource(dict(
        xs=[[r['dp_x0'], r['dp_x1']] for r in cross],
        ys=[[r['dp_y0'], r['dp_y1']] for r in cross],
        line_color=[lighten(PALETTE[r['palette_index'] % DEFAULT_COLORS]) for r in cross],
        alpha=[r['alpha'] for r in cross],
        label=[r['label'] for r in cross],
        palette_index=[r['palette_index'] for r in cross],
    ))

    # every (target chromosome, comparison chromosome) cell of the dotplot
    # grid, including pairs with zero links -- the axis rulers can only
    # select one chromosome at a time, so this is the only way to make an
    # empty cell (no synteny between that pair at all) clickable. Rendered
    # invisible (see build_page()); its only job is hit-testing for tap/hover.
    cell_xs, cell_ys, cell_q, cell_s = [], [], [], []
    for qname, qsize in ds.query_chroms:
        qx0 = ds.dp_query_offsets[qname]
        for sname, ssize in ds.subject_chroms:
            sy0 = ds.dp_subject_offsets[sname]
            cell_xs.append([qx0, qx0 + qsize, qx0 + qsize, qx0])
            cell_ys.append([sy0, sy0, sy0 + ssize, sy0 + ssize])
            cell_q.append(qname)
            cell_s.append(sname)
    cell_src = ColumnDataSource(dict(xs=cell_xs, ys=cell_ys, q_name=cell_q, s_name=cell_s))

    return (
        ColumnDataSource(dict(xs=q_xs, ys=q_ys, fill_color=q_fill, name=q_name, size_label=q_size, group=q_group)),
        ColumnDataSource(dict(xs=s_xs, ys=s_ys, fill_color=s_fill, name=s_name, size_label=s_size,
                               group=s_group, palette_index=s_idx)),
        seg_src,
        ColumnDataSource(dict(xs=grid_xs, ys=grid_ys)),
        ColumnDataSource(dict(x=q_label_x, y=q_label_y, text=q_label_text)),
        ColumnDataSource(dict(x=s_label_x, y=s_label_y, text=s_label_text)),
        cell_src,
        total_x, total_y, rx, ry,
    )


# JS port of the linear-layout math above (an angle-free version of
# circular_offsets + the bar/ribbon polygon construction), generalized so
# either side can be the "pivot" (the chromosome that was clicked, drawn
# alone on its own row) with the other side packed -- see module docstring
# for why this duplication is here instead of only in Python.
SHARED_JS = r"""
window.SYN = window.SYN || {};
// mode is 'pivot' (a ring wedge or dotplot ruler was clicked -- pivotSide/
// pivotName apply) or 'pair' (a dotplot grid cell was clicked -- pairTarget/
// pairSubject apply). The color/min-block-anchors spinners need to know
// which of SYN.buildDetailData/SYN.buildPairDetailData to re-run.
// dpQueryOffsets/dpSubjectOffsets are the dotplot's *current* per-chromosome
// axis offsets (natural order at load, replaced wholesale by the "order by
// similarity" toggle -- see SYN.applyDotplotOrder) -- every dotplot redraw
// that isn't a reorder itself (recoloring, the min-block-anchors filter)
// reads these rather than recomputing an order, so it always draws whatever
// order is currently active without needing to know which one that is.
// targetLabel/comparisonLabel are the viewer-editable display names (see
// SYN.applyLabels) -- every title/bar-label that would otherwise show the
// pipeline's literal "target"/"comparison" role tags reads from these
// instead, seeded to those same literal tags at load (see build_page()'s
// seed script) so nothing changes on screen until a viewer actually types
// something into the label inputs.
SYN.state = {
    pivotSide: null, pivotName: null, mode: null, pairTarget: null, pairSubject: null,
    dpQueryOffsets: null, dpSubjectOffsets: null,
    targetLabel: null, comparisonLabel: null,
    // the zoom panel's own x_range/y_range are rewritten on every click (see
    // SYN.applyDetail below), unlike the ring/dotplot's fixed ranges -- so
    // its double-click-to-reset target has to track whatever was most
    // recently set programmatically, not a single fixed value baked in at
    // page load (see build_page()'s three DoubleTap handlers)
    detailRange: null,
};

// The label TextInputs are free-text (see build_page()'s target_label_input/
// comparison_label_input) and, unlike chromosome names, are never validated
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

// comparison wedges are never filtered by the min-block-anchors slider (only
// ribbons are), so recoloring them only ever needs the current color count
SYN.recolorSubjectWedges = function(k, subjectSource) {
    subjectSource.data['fill_color'] = subjectSource.data['palette_index'].map((idx) => SYN.paletteColor(idx, k));
    subjectSource.change.emit();
};

// SYN.data.ribbons holds EVERY ribbon (see build_overview_sources() -- the
// embedded links are always generated at MBA_SLIDER_MIN, the least strict
// slider setting), so both filtering by min-block-anchors and recoloring by
// the current color count happen together here, client-side, on every move
// of either slider -- rebuilding from the same master list rather than
// mutating whatever happens to be currently displayed keeps the two sliders
// from stepping on each other (e.g. recoloring after a filter must recolor
// only the still-visible subset, not silently undo the filter).
SYN.buildOverviewRibbons = function(minScore, k) {
    const records = SYN.data.ribbons.filter((r) => r.score >= minScore);
    return {
        xs: records.map((r) => r.xs),
        ys: records.map((r) => r.ys),
        fill_color: records.map((r) => SYN.lighten(SYN.paletteColor(r.palette_index, k))),
        alpha: records.map((r) => r.alpha),
        label: records.map((r) => r.label),
        palette_index: records.map((r) => r.palette_index),
    };
};

SYN.applyOverviewRibbons = function(minScore, k, ribbonSource) {
    ribbonSource.data = SYN.buildOverviewRibbons(minScore, k);
    ribbonSource.change.emit();
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
            const x0 = queryOffsets[l.q_chrom] + l.q_start;
            const x1 = queryOffsets[l.q_chrom] + l.q_end;
            const yLo = subjectOffsets[l.s_chrom] + l.s_start;
            const yHi = subjectOffsets[l.s_chrom] + l.s_end;
            // '-' (inverted) blocks run the comparison span the opposite
            // direction from the target span, same convention as the
            // Python-side initial render (see build_overview_sources())
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

// The default order: target reads left-to-right in its natural (file) order
// -- same as the ring and as Dataset.__init__'s dp_query_offsets -- and
// comparison reads top-to-bottom in that SAME natural order, which means
// reversing it before laying it out bottom-to-top (see Dataset.__init__'s
// dp_subject_offsets for why: linear_offsets gives the first name the
// lowest -- bottommost -- position, so the first name has to be last in the
// list to end up at the top).
SYN.dpNaturalOrder = function() {
    return {
        queryOrder: SYN.data.queryNames.slice(),
        subjectOrder: SYN.data.subjectNames.slice().reverse(),
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
// (queryOrder, subjectOrder) pair -- a JS port of build_dotplot_sources(),
// kept in sync deliberately (see module docstring on why this geometry
// exists twice) since the order can now change after page load, which
// Python's one-time render can't react to. rx/ry/totalX/totalY never change
// under a reorder (same chromosomes, same sizes, just shuffled), so they
// come from the initial Python render rather than being recomputed here.
SYN.buildDotplotLayout = function(queryOrder, subjectOrder, k) {
    const rx = SYN.data.dpRx, ry = SYN.data.dpRy;
    const totalX = SYN.data.dpTotalX, totalY = SYN.data.dpTotalY;
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
        sName.push(name); sSize.push((size / 1e6).toFixed(2) + ' Mb'); sGroup.push('comparison');
        sIdx.push(SYN.data.subjectColorIndex[name]);
    }

    const gridXs = [], gridYs = [];
    for (const name of queryOrder) {
        const x0 = queryOffsets[name];
        gridXs.push([x0, x0]); gridYs.push([-ry, totalY]);
    }
    const lastQ = queryOrder[queryOrder.length - 1];
    gridXs.push([queryOffsets[lastQ] + SYN.data.querySizes[lastQ], queryOffsets[lastQ] + SYN.data.querySizes[lastQ]]);
    gridYs.push([-ry, totalY]);
    for (const name of subjectOrder) {
        const y0 = subjectOffsets[name];
        gridXs.push([-rx, totalX]); gridYs.push([y0, y0]);
    }
    const lastS = subjectOrder[subjectOrder.length - 1];
    gridXs.push([-rx, totalX]);
    gridYs.push([subjectOffsets[lastS] + SYN.data.subjectSizes[lastS], subjectOffsets[lastS] + SYN.data.subjectSizes[lastS]]);

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

    const cellXs = [], cellYs = [], cellQ = [], cellS = [];
    for (const qn of queryOrder) {
        const qx0 = queryOffsets[qn], qsize = SYN.data.querySizes[qn];
        for (const sn of subjectOrder) {
            const sy0 = subjectOffsets[sn], ssize = SYN.data.subjectSizes[sn];
            cellXs.push([qx0, qx0 + qsize, qx0 + qsize, qx0]);
            cellYs.push([sy0, sy0, sy0 + ssize, sy0 + ssize]);
            cellQ.push(qn); cellS.push(sn);
        }
    }

    return {
        queryOffsets, subjectOffsets,
        queryRuler: {xs: qXs, ys: qYs, fill_color: queryOrder.map(() => SYN.data.targetGrey),
                     name: qName, size_label: qSize, group: qGroup},
        subjectRuler: {xs: sXs, ys: sYs, fill_color: sIdx.map((idx) => SYN.paletteColor(idx, k)),
                       name: sName, size_label: sSize, group: sGroup, palette_index: sIdx},
        grid: {xs: gridXs, ys: gridYs},
        queryLabels: {x: qLabelX, y: qLabelY, text: qLabelText},
        subjectLabels: {x: sLabelX, y: sLabelY, text: sLabelText},
        cells: {xs: cellXs, ys: cellYs, q_name: cellQ, s_name: cellS},
    };
};

// Applies a full reorder: rebuilds every order-dependent dotplot source (via
// buildDotplotLayout) and updates SYN.state's offsets so every later
// recolor/refilter (color spinner, palette dropdown, min-block-anchors) draws
// against the new order without needing to know a reorder happened.
SYN.applyDotplotOrder = function(queryOrder, subjectOrder, k, minScore, sources) {
    const layout = SYN.buildDotplotLayout(queryOrder, subjectOrder, k);
    SYN.state.dpQueryOffsets = layout.queryOffsets;
    SYN.state.dpSubjectOffsets = layout.subjectOffsets;
    sources.query.data = layout.queryRuler; sources.query.change.emit();
    sources.subject.data = layout.subjectRuler; sources.subject.change.emit();
    sources.grid.data = layout.grid; sources.grid.change.emit();
    sources.queryLabel.data = layout.queryLabels; sources.queryLabel.change.emit();
    sources.subjectLabel.data = layout.subjectLabels; sources.subjectLabel.change.emit();
    sources.cell.data = layout.cells; sources.cell.change.emit();
    SYN.applyDotplotSegmentsForCurrentLayout(minScore, k, sources.segment);
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

SYN.emptyDetail = function(title) {
    return {
        bars: {xs: [], ys: [], fill_color: []},
        ribbons: {xs: [], ys: [], fill_color: [], alpha: [], label: []},
        labels: {x: [], y: [], text: [], color: []},
        x_range: [-1, 1], y_range: [-0.5, 1.5],
        title: title || 'Click any chromosome wedge on the left to zoom in.',
    };
};

// pivotSide is 'subject' (comparison) or 'query' (target) -- whichever side
// was clicked. That chromosome is drawn alone on its own row; the other
// side's chromosomes involved in its links are packed on the opposite row.
// Comparison-associated bars/ribbons always end up on the top row, target-
// associated ones on the bottom row, regardless of which side is the pivot.
SYN.buildDetailData = function(pivotSide, pivotName, k, minScore) {
    if (!pivotName) { return SYN.emptyDetail(); }
    const isPivotSubject = pivotSide === 'subject';
    const links = ((isPivotSubject ? SYN.data.linksByComparison[pivotName] : SYN.data.linksByQuery[pivotName]) || [])
        .filter((l) => l.score >= minScore);
    const pivotSize = isPivotSubject ? SYN.data.subjectSizes[pivotName] : SYN.data.querySizes[pivotName];
    const pivotColor = isPivotSubject ? SYN.paletteColor(SYN.data.subjectColorIndex[pivotName], k) : SYN.data.targetGrey;
    const pivotKind = isPivotSubject ? SYN.state.comparisonLabel : SYN.state.targetLabel;
    const otherKind = isPivotSubject ? SYN.state.targetLabel : SYN.state.comparisonLabel;

    const barH = 0.32, topY = 1.2, botY = 0.0;
    const pivotRowY = isPivotSubject ? topY : botY;
    const otherRowY = isPivotSubject ? botY : topY;

    if (links.length === 0) {
        const d = SYN.emptyDetail(`${pivotName}: no links found`);
        d.bars = {xs: [[0, pivotSize, pivotSize, 0]],
                  ys: [[pivotRowY, pivotRowY, pivotRowY + barH, pivotRowY + barH]],
                  fill_color: [pivotColor]};
        d.labels = {x: [pivotSize / 2], y: [pivotRowY + barH / 2], text: [pivotName],
                    color: [SYN.textColorFor(pivotColor)]};
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
    }
    barXs.push([0, pivotSize, pivotSize, 0]);
    barYs.push([pivotRowY, pivotRowY, pivotRowY + barH, pivotRowY + barH]);
    barFill.push(pivotColor);
    labelX.push(pivotSize / 2);
    labelY.push(pivotRowY + barH / 2);
    labelText.push(pivotName);
    labelColor.push(SYN.textColorFor(pivotColor));

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
        x_range: [-xMax * 0.03, xMax * 1.05],
        y_range: [botY - 0.35, topY + barH + 0.35],
        title: `${pivotName} (${pivotKind}) — ${links.length} link(s) vs ${otherNames.length} ${otherKind} chromosome(s)`,
    };
};

// One dotplot grid cell (a single target chromosome x single comparison
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

    if (links.length === 0) {
        return {
            bars: {xs: barXs, ys: barYs, fill_color: barFill},
            ribbons: {xs: [], ys: [], fill_color: [], alpha: [], label: []},
            labels: {x: labelX, y: labelY, text: labelText, color: labelColor},
            x_range: baseX, y_range: baseY,
            title: `${targetName} (${SYN.state.targetLabel}) × ${subjectName} (${SYN.state.comparisonLabel}) -- no links`,
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
        x_range: baseX, y_range: baseY,
        title: `${targetName} (${SYN.state.targetLabel}) × ${subjectName} (${SYN.state.comparisonLabel}) -- ${links.length} link(s)`,
    };
};

// mirrors format_link_label() in this file's Python -- mean_identity/
// anchor_density are null for links.tsv files predating those columns, so
// that line is simply omitted rather than printing "null"
SYN.formatLinkLabel = function(l) {
    let stats = `${l.score} protein alignment(s)  ·  orientation=${l.orientation}`;
    if (l.mean_identity !== null && l.mean_identity !== undefined
        && l.anchor_density !== null && l.anchor_density !== undefined) {
        stats += `<br>avg identity: ${(l.mean_identity * 100).toFixed(1)}%`
            + `  ·  density: ${l.anchor_density.toFixed(1)} anchors/Mb`;
    }
    return `${l.q_chrom}:${l.q_start.toLocaleString()}-${l.q_end.toLocaleString()}`
        + ` ↔ ${l.s_chrom}:${l.s_start.toLocaleString()}-${l.s_end.toLocaleString()}`
        + `<br>${stats}`;
};

SYN.applyDetail = function(d, barSource, ribbonSource, labelSource, detailFig) {
    barSource.data = d.bars; barSource.change.emit();
    ribbonSource.data = d.ribbons; ribbonSource.change.emit();
    labelSource.data = d.labels; labelSource.change.emit();
    detailFig.x_range.start = d.x_range[0]; detailFig.x_range.end = d.x_range[1];
    detailFig.y_range.start = d.y_range[0]; detailFig.y_range.end = d.y_range[1];
    detailFig.title.text = d.title;
    SYN.state.detailRange = {x0: d.x_range[0], x1: d.x_range[1], y0: d.y_range[0], y1: d.y_range[1]};
};

// Re-renders whatever is currently shown in the zoom panel (a pivot, a
// pair, or nothing) at a new color count / min-block-anchors threshold --
// shared by the color-count spinner, the palette dropdown, and the min-
// block-anchors spinner (see build_page()), so which of buildDetailData/
// buildPairDetailData applies doesn't need to be duplicated in each of
// their three CustomJS callbacks.
SYN.refreshDetail = function(k, minScore, barSource, ribbonSource, labelSource, detailFig) {
    if (SYN.state.mode === 'pair') {
        const d = SYN.buildPairDetailData(SYN.state.pairTarget, SYN.state.pairSubject, k, minScore);
        SYN.applyDetail(d, barSource, ribbonSource, labelSource, detailFig);
    } else if (SYN.state.pivotName) {
        const d = SYN.buildDetailData(SYN.state.pivotSide, SYN.state.pivotName, k, minScore);
        SYN.applyDetail(d, barSource, ribbonSource, labelSource, detailFig);
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

// Exports a figure as a real, standalone .svg file -- no PDF conversion
// library needed at all, since output_backend="svg" already gives every
// figure a genuine vector <svg> DOM tree (see build_page()) rather than only
// ever being able to rasterize a canvas. A user who actually wants a PDF (or
// PNG, or anything else) can convert this in Inkscape, Illustrator, or any
// other vector tool -- exporting the SVG directly here keeps this page free
// of any runtime dependency beyond Bokeh itself.
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
    const blob = new Blob([svgText], {type: 'image/svg+xml;charset=utf-8'});
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
};

// Builds a save button's download filename from whatever labels are
// currently active, instead of a fixed generic name, so an exported file is
// self-describing without the viewer having to rename it -- sanitized since
// the labels are free-text input and filenames can't contain arbitrary
// characters (path separators in particular).
SYN.buildExportFilename = function(panel) {
    const sanitize = (s) => (s || '').replace(/[^a-zA-Z0-9_-]+/g, '_').replace(/^_+|_+$/g, '');
    const t = sanitize(SYN.state.targetLabel) || 'target';
    const r = sanitize(SYN.state.comparisonLabel) || 'comparison';
    return `${t}_vs_${r}_${panel}.svg`;
};

// Every title/bar-label a viewer might export reads SYN.state.targetLabel/
// comparisonLabel rather than the pipeline's literal "target"/"comparison"
// role tags (see SYN.buildDetailData/SYN.buildPairDetailData above) -- this
// is what actually changes those two values and refreshes the handful of
// sources/titles that don't go through a rebuild-on-every-redraw path the
// way the zoom panel does (see build_page()'s label-input callback, which
// also calls SYN.refreshDetail and SYN.applyStats after this). The 'group'
// field feeds each ruler/wedge source's hover tooltip ("Group": "@group")
// only -- nothing else reads it, so overwriting every entry is safe.
SYN.applyLabels = function(targetLabel, comparisonLabel, ctx) {
    SYN.state.targetLabel = targetLabel;
    SYN.state.comparisonLabel = comparisonLabel;
    ctx.overview.title.text = `${targetLabel} (target) vs ${comparisonLabel} (comparison)`;
    ctx.dotplotFig.title.text = `${targetLabel} vs ${comparisonLabel} -- click a band or a grid square to zoom`;
    if (ctx.headerDiv) {
        ctx.headerDiv.text = `<h2>${SYN.escapeHtml(targetLabel)} vs ${SYN.escapeHtml(comparisonLabel)} `
            + `-- interactive synteny</h2>`;
    }
    for (const src of ctx.queryLikeSources) {
        src.data.group = src.data.group.map(() => targetLabel);
        src.change.emit();
    }
    for (const src of ctx.subjectLikeSources) {
        src.data.group = src.data.group.map(() => comparisonLabel);
        src.change.emit();
    }
};

// The interactive HTML's replacement for the (removed) static plots' stats
// inset -- see compute_alignment_stats.py and build_page()'s --stats. Lives
// as a normal page element rather than baked into any one exportable figure,
// since export is per-panel now (ring, zoom, or dotplot), not "the whole
// image" the way a static render was.
SYN.buildStatsHtml = function(targetLabel, comparisonLabel) {
    const s = SYN.data.alignmentStats;
    if (!s) { return ''; }
    const pct = (n) => s.proteome_total ? (n / s.proteome_total * 100).toFixed(1) : '0.0';
    const row = (label, aligned, identity) => `
        <tr><td>${SYN.escapeHtml(label)}</td><td>${aligned.toLocaleString()} (${pct(aligned)}%)</td>
            <td>${(identity * 100).toFixed(1)}%</td></tr>`;
    return `
        <div style="font-size:13px;color:#333;border:1px solid #ddd;border-radius:6px;
                    padding:8px 12px;display:inline-block;margin:4px 0 8px 0">
            <b>Alignment summary</b> -- ${s.proteome_total.toLocaleString()} protein(s) in the input proteome
            <table style="border-collapse:collapse;margin-top:4px">
                <tr style="color:#666"><th style="text-align:left;padding-right:12px"></th>
                    <th style="text-align:right;padding-right:12px">aligned</th>
                    <th style="text-align:right">avg identity</th></tr>
                ${row(targetLabel, s.query_aligned, s.query_mean_identity)}
                ${row(comparisonLabel, s.subject_aligned, s.subject_mean_identity)}
            </table>
        </div>`;
};

SYN.applyStats = function(targetLabel, comparisonLabel, statsDiv) {
    statsDiv.text = SYN.buildStatsHtml(targetLabel, comparisonLabel);
};

"""


def build_page(ds, query_name, subject_name, query_subtitle=None, subject_subtitle=None,
               alignment_stats=None):
    query_names = [n for n, _ in ds.query_chroms]
    subject_names = [n for n, _ in ds.subject_chroms]
    subject_index = {name: i for i, (name, _) in enumerate(ds.subject_chroms)}
    query_index = {name: i for i, (name, _) in enumerate(ds.query_chroms)}

    q_src, s_src, r_src, label_src, ribbon_records = build_overview_sources(ds, subject_index, query_index)

    lim = OUTER_R + 0.22
    # no 'tap' here -- that shorthand adds its own generic, unrestricted
    # TapTool, which would sit in the toolbar (visible) alongside the
    # renderer-scoped, hidden one added below via add_tools(), showing two
    # redundant Tap buttons for the same gesture. The renderer-scoped one
    # is the one actually driving click-to-zoom (see .selected.js_on_change
    # wiring below) and is the only one that needs to exist at all.
    overview = figure(width=620, height=620, match_aspect=True,
                       x_range=Range1d(-lim, lim), y_range=Range1d(-lim, lim),
                       title=f"{query_name} (target) vs {subject_name} (comparison)",
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

    # follow_mouse -- see detail_fig's identical HoverTool below for why: a
    # ribbon can span most of the ring, so the default snap-to-center
    # tooltip can land outside the current (zoomed/panned) view
    overview.add_tools(HoverTool(renderers=[ribbon_renderer], tooltips="@label{safe}",
                                  point_policy='follow_mouse'))
    overview.add_tools(HoverTool(renderers=[query_renderer, subject_renderer],
                                  tooltips=[("Chromosome", "@name"), ("Size", "@size_label"),
                                            ("Group", "@group")]))
    # both sides are clickable: either wedge triggers the same zoom behavior.
    # visible=False hides its toolbar button -- clicking a wedge doesn't
    # need one, it's not a mode you toggle on/off -- but a tool added via
    # add_tools() (unlike one named in the tools= string above) is never
    # auto-activated, hidden or not, so active_tap has to be set explicitly
    # or the hidden tool would sit there never actually receiving clicks
    overview_tap = TapTool(renderers=[subject_renderer, query_renderer], visible=False)
    overview.add_tools(overview_tap)
    overview.toolbar.active_tap = overview_tap
    # double-click to reset pan/zoom -- a shortcut for the toolbar's own
    # Reset button, since this ring never reprograms its own range after
    # creation, the reset target is just its own original (fixed) extent
    overview.js_on_event(DoubleTap, CustomJS(args=dict(fig=overview), code=f"""
        fig.x_range.start = {-lim}; fig.x_range.end = {lim};
        fig.y_range.start = {-lim}; fig.y_range.end = {lim};
    """))

    # Whole-genome dotplot: target along x, comparison along y, each block
    # drawn as the line segment its own (query, subject) span traces out --
    # the standard synteny-dotplot convention (a diagonal streak per
    # collinear block, anti-diagonal for an inversion). The chromosome-ruler
    # strips along each axis are clickable, wired into the exact same
    # zoom-pivot mechanism the ring wedges use below (same 'name' field, same
    # tap callback) -- from SYN.buildDetailData's perspective, a click is a
    # click regardless of which panel it came from.
    (dp_q_src, dp_s_src, dp_seg_src, dp_grid_src,
     dp_q_label_src, dp_s_label_src, dp_cell_src, dp_total_x, dp_total_y, dp_rx, dp_ry) = \
        build_dotplot_sources(ds, subject_index, ribbon_records)

    # *3 (not just enough room for the ruler strip itself) leaves space for
    # the chromosome-name labels drawn just outside each ruler (see
    # build_dotplot_sources()'s q_label/s_label -- text_align='right' for the
    # comparison labels means they extend further left from their anchor, so
    # the range has to clear that too, not just the ruler). Confirmed by
    # rendering against real data: *2.2 and *2.6 clip a 5-character name like
    # "chr13" down to "13"/"hr13"; *3 leaves every label fully visible with a
    # little room to spare while still keeping the dead margin well under
    # half of what a much larger multiplier would waste.
    # no 'tap' here -- see overview's identical figure(...) above for why
    dotplot_fig = figure(width=620, height=620,
                          x_range=Range1d(-dp_rx * 3, dp_total_x * 1.02),
                          y_range=Range1d(-dp_ry * 3, dp_total_y * 1.02),
                          title="Whole-genome dotplot -- click a band or a grid square to zoom",
                          tools="pan,wheel_zoom,reset",
                          output_backend="svg")
    dotplot_fig.axis.visible = False
    dotplot_fig.grid.visible = False
    dotplot_fig.toolbar.logo = None

    # invisible hit-target layer, one patch per (target chromosome, comparison
    # chromosome) grid cell -- drawn first (bottom of the stack) so it never
    # visually covers the gridlines/segments/rulers above it. This is what
    # makes every cell clickable for the pairwise zoom below, including
    # cells with zero links, which the axis rulers alone can't reach.
    dp_cell_renderer = dotplot_fig.patches('xs', 'ys', source=dp_cell_src, fill_color=None,
                                            fill_alpha=0, line_color=None)
    dotplot_fig.multi_line('xs', 'ys', source=dp_grid_src, line_color='#DDDDDD', line_width=1)
    dp_segment_renderer = dotplot_fig.multi_line('xs', 'ys', source=dp_seg_src, line_color='line_color',
                                                  line_alpha='alpha', line_width=2)
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
    # one TapTool covering every clickable renderer here (rather than one
    # per renderer group, as an earlier version had) -- a single tool can
    # span renderers backed by different data sources just fine, each
    # click still only touches whichever renderer/glyph was actually hit,
    # and it means only one tool's active_tap needs setting below.
    # visible=False -- see overview's identical TapTool above for why;
    # grid cells get no hover tooltip (see SYN.buildPairDetailData below):
    # the axis rulers' own hover already names the chromosome under the
    # cursor and a per-cell tooltip on top of it was judged redundant/noisy
    # given how densely the cells tile the plot
    dotplot_tap = TapTool(renderers=[dp_query_renderer, dp_subject_renderer, dp_cell_renderer], visible=False)
    dotplot_fig.add_tools(dotplot_tap)
    dotplot_fig.toolbar.active_tap = dotplot_tap
    # double-click to reset -- see overview's identical handler above; the
    # dotplot's range is likewise fixed at creation and never reprogrammed
    # (reordering chromosomes only changes their offsets within that same
    # fixed total span, see SYN.buildDotplotLayout)
    dotplot_fig.js_on_event(DoubleTap, CustomJS(args=dict(fig=dotplot_fig), code=f"""
        fig.x_range.start = {-dp_rx * 3}; fig.x_range.end = {dp_total_x * 1.02};
        fig.y_range.start = {-dp_ry * 3}; fig.y_range.end = {dp_total_y * 1.02};
    """))

    detail_bar_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[]))
    detail_rib_src = ColumnDataSource(dict(xs=[], ys=[], fill_color=[], alpha=[], label=[]))
    detail_label_src = ColumnDataSource(dict(x=[], y=[], text=[], color=[]))
    detail_fig = figure(width=446, height=302, x_axis_label='position (Mb)',
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
    detail_fig.text('x', 'y', source=detail_label_src, text_align='center', text_baseline='middle',
                     text_font_size='9px', text_color='color')
    # point_policy='follow_mouse' (default is 'snap_to_data', anchored to
    # the ribbon's geometric center) -- a ribbon can span the whole panel,
    # so anchoring to its center means the tooltip can render outside the
    # current view once someone's zoomed into one end of it
    detail_fig.add_tools(HoverTool(renderers=[detail_ribbon_renderer], tooltips="@label{safe}",
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

    # default to the literal role tags (today's look, unchanged) -- a viewer
    # can retype either to the actual species/genome name before exporting a
    # panel, see SYN.applyLabels
    target_label_input = TextInput(title="Target label", value=query_name, width=220)
    comparison_label_input = TextInput(title="Comparison label", value=subject_name, width=220)

    palette_select = Select(title="Color palette", value=DEFAULT_PALETTE_NAME,
                             options=list(PALETTES.keys()), width=220)
    color_spinner = Spinner(title="Comparison colors (discrete, cycling)", low=1, high=MAX_COLORS,
                             step=1, value=DEFAULT_COLORS, width=220)
    # high is the actual largest block size in this dataset, not a fixed
    # constant -- MBA_SLIDER_MIN is still a real floor (the embedded links
    # were generated with --min_block_anchors that low, see
    # build_synteny.nf's *_FOR_SLIDER processes, so nothing smaller was ever
    # computed to reveal), but there's no reason to cap the top end at some
    # arbitrary round number when a Spinner lets someone type any value
    max_observed_score = max((r['score'] for r in ribbon_records), default=MBA_SLIDER_MIN)
    mba_spinner = Spinner(title="Min block anchors", low=MBA_SLIDER_MIN, high=max_observed_score,
                           step=1, value=MBA_SLIDER_MIN, width=220)
    reset_btn = Button(label="✕ clear zoom", button_type="default", width=140)
    save_ring_btn = Button(label="⬇ save ring (SVG)", button_type="default", width=160)
    save_zoom_btn = Button(label="⬇ save zoom (SVG)", button_type="default", width=160)
    save_dotplot_btn = Button(label="⬇ save dotplot (SVG)", button_type="default", width=160)
    # off by default: natural (file/karyotype) order is what most users
    # recognize their chromosomes by, and reordering only pays off when the
    # two genomes are close enough that a near-1:1 correspondence exists to
    # reveal in the first place (see SYN.computeSimilarityOrder's docstring)
    order_toggle = Toggle(label="Order chromosomes by similarity", active=False,
                           button_type="default", width=160)

    save_ring_btn.js_on_click(CustomJS(args=dict(fig=overview), code="""
        SYN.exportFigureAsSVG(fig, SYN.buildExportFilename('ring'));
    """))
    save_zoom_btn.js_on_click(CustomJS(args=dict(fig=detail_fig), code="""
        SYN.exportFigureAsSVG(fig, SYN.buildExportFilename('zoom'));
    """))
    save_dotplot_btn.js_on_click(CustomJS(args=dict(fig=dotplot_fig), code="""
        SYN.exportFigureAsSVG(fig, SYN.buildExportFilename('dotplot'));
    """))

    # Rebuilds every order-dependent dotplot source for whichever order is
    # now active -- SYN.applyDotplotOrder also updates SYN.state's offsets,
    # so the color/palette/min-block-anchors callbacks below (which redraw
    # the dotplot too, via SYN.applyDotplotSegmentsForCurrentLayout) pick up
    # the new order automatically without needing to know a reorder happened.
    order_toggle_callback = CustomJS(args=dict(
        color_spinner=color_spinner, mba_spinner=mba_spinner,
        dp_query_source=dp_q_src, dp_subject_source=dp_s_src, dp_grid_source=dp_grid_src,
        dp_q_label_source=dp_q_label_src, dp_s_label_source=dp_s_label_src,
        dp_cell_source=dp_cell_src, dp_segment_source=dp_seg_src,
    ), code="""
        const order = cb_obj.active ? SYN.computeSimilarityOrder() : SYN.dpNaturalOrder();
        SYN.applyDotplotOrder(order.queryOrder, order.subjectOrder, color_spinner.value, mba_spinner.value, {
            query: dp_query_source, subject: dp_subject_source, grid: dp_grid_source,
            queryLabel: dp_q_label_source, subjectLabel: dp_s_label_source,
            cell: dp_cell_source, segment: dp_segment_source,
        });
    """)
    order_toggle.js_on_click(order_toggle_callback)

    # Tapping a chromosome region -- a ring wedge, a dotplot ruler cell, or a
    # dotplot grid cell -- clears every OTHER clickable source's selection,
    # so exactly one chromosome (or pair) is ever "active" across all three
    # panels at once. other_sources is a list because there are now five
    # mutually-exclusive click sources (ring x2, dotplot ruler x2, dotplot
    # grid), not two -- each needs to clear the other four. Two side groups
    # ('subject'/'query') rather than one callback per source because
    # cb_obj here is the Selection model that changed, not the
    # ColumnDataSource itself, so there's nothing reliable to branch on
    # inside a single shared callback.
    def make_tap_callback(this_source, other_sources, side):
        return CustomJS(args=dict(
            this_source=this_source, other_sources=other_sources, side=side,
            color_spinner=color_spinner, mba_spinner=mba_spinner, bar_source=detail_bar_src,
            ribbon_source=detail_rib_src, label_source=detail_label_src, detail_fig=detail_fig,
        ), code="""
            const idx = this_source.selected.indices;
            if (idx.length === 0) { return; }
            for (const other of other_sources) { other.selected.indices = []; }
            const name = this_source.data['name'][idx[idx.length - 1]];
            SYN.state.mode = 'pivot';
            SYN.state.pivotSide = side;
            SYN.state.pivotName = name;
            const d = SYN.buildDetailData(side, name, color_spinner.value, mba_spinner.value);
            SYN.applyDetail(d, bar_source, ribbon_source, label_source, detail_fig);
        """)

    s_src.selected.js_on_change(
        'indices', make_tap_callback(s_src, [q_src, dp_s_src, dp_q_src, dp_cell_src], 'subject'))
    q_src.selected.js_on_change(
        'indices', make_tap_callback(q_src, [s_src, dp_s_src, dp_q_src, dp_cell_src], 'query'))
    dp_s_src.selected.js_on_change(
        'indices', make_tap_callback(dp_s_src, [s_src, q_src, dp_q_src, dp_cell_src], 'subject'))
    dp_q_src.selected.js_on_change(
        'indices', make_tap_callback(dp_q_src, [s_src, q_src, dp_s_src, dp_cell_src], 'query'))

    # Tapping a dotplot grid cell zooms the detail panel into exactly that
    # (target chromosome, comparison chromosome) pair -- unlike the pivot
    # callback above, this always shows exactly two bars, and handles zero
    # links as a normal result (see SYN.buildPairDetailData) rather than
    # refusing the click, since "no synteny here" is itself the answer for
    # an empty cell.
    pair_tap_callback = CustomJS(args=dict(
        this_source=dp_cell_src, other_sources=[s_src, q_src, dp_s_src, dp_q_src],
        color_spinner=color_spinner, mba_spinner=mba_spinner, bar_source=detail_bar_src,
        ribbon_source=detail_rib_src, label_source=detail_label_src, detail_fig=detail_fig,
    ), code="""
        const idx = this_source.selected.indices;
        if (idx.length === 0) { return; }
        for (const other of other_sources) { other.selected.indices = []; }
        const i = idx[idx.length - 1];
        const targetName = this_source.data['q_name'][i];
        const subjectName = this_source.data['s_name'][i];
        SYN.state.mode = 'pair';
        SYN.state.pairTarget = targetName;
        SYN.state.pairSubject = subjectName;
        const d = SYN.buildPairDetailData(targetName, subjectName, color_spinner.value, mba_spinner.value);
        SYN.applyDetail(d, bar_source, ribbon_source, label_source, detail_fig);
    """)
    dp_cell_src.selected.js_on_change('indices', pair_tap_callback)

    color_callback = CustomJS(args=dict(
        subject_source=s_src, dp_subject_source=dp_s_src,
        overview_ribbon_source=r_src, dotplot_segment_source=dp_seg_src, mba_spinner=mba_spinner,
        bar_source=detail_bar_src, ribbon_source=detail_rib_src,
        label_source=detail_label_src, detail_fig=detail_fig,
    ), code="""
        const k = cb_obj.value;
        SYN.recolorSubjectWedges(k, subject_source);
        SYN.recolorSubjectWedges(k, dp_subject_source);
        SYN.applyOverviewRibbons(mba_spinner.value, k, overview_ribbon_source);
        SYN.applyDotplotSegmentsForCurrentLayout(mba_spinner.value, k, dotplot_segment_source);
        SYN.refreshDetail(k, mba_spinner.value, bar_source, ribbon_source, label_source, detail_fig);
    """)
    color_spinner.js_on_change('value', color_callback)

    mba_callback = CustomJS(args=dict(
        overview_ribbon_source=r_src, dotplot_segment_source=dp_seg_src, color_spinner=color_spinner,
        bar_source=detail_bar_src, ribbon_source=detail_rib_src,
        label_source=detail_label_src, detail_fig=detail_fig,
    ), code="""
        const minScore = cb_obj.value;
        const k = color_spinner.value;
        SYN.applyOverviewRibbons(minScore, k, overview_ribbon_source);
        SYN.applyDotplotSegmentsForCurrentLayout(minScore, k, dotplot_segment_source);
        SYN.refreshDetail(k, minScore, bar_source, ribbon_source, label_source, detail_fig);
    """)
    mba_spinner.js_on_change('value', mba_callback)

    # Swaps the active color array itself (SYN.data.palette -- every recolor
    # function reads from it by reference, so nothing else needs to change)
    # rather than the number of colors cycled through, which is what
    # color_spinner controls -- the two are independent and compose (e.g.
    # "Pastel" at a count of 4).
    palette_callback = CustomJS(args=dict(
        subject_source=s_src, dp_subject_source=dp_s_src,
        overview_ribbon_source=r_src, dotplot_segment_source=dp_seg_src,
        color_spinner=color_spinner, mba_spinner=mba_spinner,
        bar_source=detail_bar_src, ribbon_source=detail_rib_src,
        label_source=detail_label_src, detail_fig=detail_fig,
    ), code="""
        SYN.data.palette = SYN.data.palettes[cb_obj.value];
        const k = color_spinner.value;
        SYN.recolorSubjectWedges(k, subject_source);
        SYN.recolorSubjectWedges(k, dp_subject_source);
        SYN.applyOverviewRibbons(mba_spinner.value, k, overview_ribbon_source);
        SYN.applyDotplotSegmentsForCurrentLayout(mba_spinner.value, k, dotplot_segment_source);
        SYN.refreshDetail(k, mba_spinner.value, bar_source, ribbon_source, label_source, detail_fig);
    """)
    palette_select.js_on_change('value', palette_callback)

    reset_callback = CustomJS(args=dict(
        subject_source=s_src, query_source=q_src, dp_subject_source=dp_s_src, dp_query_source=dp_q_src,
        dp_cell_source=dp_cell_src,
        bar_source=detail_bar_src, ribbon_source=detail_rib_src,
        label_source=detail_label_src, detail_fig=detail_fig,
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
        dp_cell_source.selected.indices = [];
        SYN.applyDetail(SYN.emptyDetail(), bar_source, ribbon_source, label_source, detail_fig);
    """)
    reset_btn.js_on_click(reset_callback)

    hint = Div(text="<p style='color:#666;font-size:13px'>Click any chromosome wedge on the "
                     "ring, or any chromosome band on the dotplot's axes, to zoom into its "
                     "links against the other genome -- both trigger the same zoom, and a "
                     "click on any of them clears the others' selection. You can also click a "
                     "single square in the dotplot grid to zoom straight into that one "
                     "chromosome pair, including squares with no alignments at all. Target "
                     "label/Comparison label relabel the genomes everywhere a title or bar shows "
                     "them, so an exported panel can show the actual species/genome name "
                     "instead of \"target\"/\"comparison\". Color palette switches the set of "
                     "colors comparison chromosomes cycle through; Comparison colors changes how "
                     "many discrete colors from it they cycle through. Min block anchors filters "
                     "out synteny blocks with fewer than that many supporting protein alignments "
                     "-- raise it to cut noise, lower it to see more (shorter, less certain) "
                     "blocks. Order chromosomes by similarity reorders both dotplot axes so "
                     "shared synteny lines up into a diagonal -- most useful when the two "
                     "genomes are closely related with a roughly 1:1 chromosome correspondence; "
                     "off by default so chromosomes stay in their familiar order. Hover any "
                     "wedge, band, or ribbon for details, and use each panel's save button to "
                     "export it as an SVG (open it in Inkscape, Illustrator, or similar to "
                     "convert to PDF/PNG/etc.).</p>")

    left_col = column(overview, save_ring_btn)
    mid_col = column(detail_fig, row(reset_btn, save_zoom_btn))
    right_col = column(dotplot_fig, row(save_dotplot_btn, order_toggle))

    # reactive (see SYN.applyLabels) -- shows the current target_label_input/
    # comparison_label_input values, not the literal "target"/"comparison"
    # role tags, once a viewer changes either
    header_title_div = Div(text=f"<h2>{html.escape(query_name)} vs {html.escape(subject_name)} "
                                 "-- interactive synteny</h2>")

    subtitle_div = None
    if query_subtitle or subject_subtitle:
        # the input file name/accession behind "target"/"comparison", so a
        # viewer doesn't have to guess which physical genome each role is --
        # this deliberately never changes when the labels above do, since
        # its whole point is mapping the fixed role to the actual input file
        subtitle_div = Div(text=(
            "<p style='color:#666;font-size:13px;margin:2px 0 8px 0'>"
            f"{html.escape(query_name)}: {html.escape(query_subtitle or '?')}<br>"
            f"{html.escape(subject_name)}: {html.escape(subject_subtitle or '?')}"
            "</p>"
        ))

    # initial render uses the default (literal role-tag) labels, matching
    # the TextInputs' own defaults -- reactive after that (see SYN.applyStats,
    # wired below), same split as header_title_div/subtitle_div above.
    # Empty/invisible if no --stats was given.
    stats_div = Div(text=format_stats_html(query_name, subject_name, alignment_stats))

    label_callback = CustomJS(args=dict(
        target_label_input=target_label_input, comparison_label_input=comparison_label_input,
        overview=overview, dotplot_fig=dotplot_fig, header_div=header_title_div, stats_div=stats_div,
        q_src=q_src, dp_q_src=dp_q_src, s_src=s_src, dp_s_src=dp_s_src,
        color_spinner=color_spinner, mba_spinner=mba_spinner,
        bar_source=detail_bar_src, ribbon_source=detail_rib_src,
        label_source=detail_label_src, detail_fig=detail_fig,
    ), code="""
        const targetLabel = target_label_input.value.trim() || 'target';
        const comparisonLabel = comparison_label_input.value.trim() || 'comparison';
        SYN.applyLabels(targetLabel, comparisonLabel, {
            overview, dotplotFig: dotplot_fig, headerDiv: header_div,
            queryLikeSources: [q_src, dp_q_src], subjectLikeSources: [s_src, dp_s_src],
        });
        SYN.refreshDetail(color_spinner.value, mba_spinner.value, bar_source, ribbon_source, label_source, detail_fig);
        SYN.applyStats(targetLabel, comparisonLabel, stats_div);
    """)
    target_label_input.js_on_change('value', label_callback)
    comparison_label_input.js_on_change('value', label_callback)

    header_children = [header_title_div]
    if subtitle_div is not None:
        header_children.append(subtitle_div)

    layout = column(
        *header_children,
        stats_div,
        row(target_label_input, comparison_label_input),
        hint,
        row(palette_select, color_spinner, mba_spinner),
        row(left_col, mid_col, right_col),
    )

    page_html = file_html(layout, CDN, title=f"{query_name} vs {subject_name} -- interactive synteny")

    def embed_link(l):
        return {'q_chrom': l['q_chrom'], 'q_start': l['q_start'], 'q_end': l['q_end'],
                's_chrom': l['s_chrom'], 's_start': l['s_start'], 's_end': l['s_end'],
                'score': l['score'], 'orientation': l['orientation'],
                'mean_identity': l.get('mean_identity'), 'anchor_density': l.get('anchor_density')}

    links_by_comparison, links_by_query = {}, {}
    for name, _ in ds.subject_chroms:
        links_by_comparison[name] = [embed_link(l) for l in ds.links if l['s_chrom'] == name]
    for name, _ in ds.query_chroms:
        links_by_query[name] = [embed_link(l) for l in ds.links if l['q_chrom'] == name]
    syn_data = {
        'palette': PALETTES[DEFAULT_PALETTE_NAME],
        'palettes': PALETTES,
        'targetGrey': TARGET_GREY,
        'targetName': query_name,
        'comparisonName': subject_name,
        'queryNames': query_names,
        'subjectNames': subject_names,
        'querySizes': ds.query_sizes,
        'subjectSizes': ds.subject_sizes,
        'subjectColorIndex': subject_index,
        'linksByComparison': links_by_comparison,
        'linksByQuery': links_by_query,
        'ribbons': ribbon_records,
        'maxScore': ds.max_score(),
        # rx/ry/totalX/totalY never change when the dotplot is reordered --
        # same chromosomes, same sizes, just shuffled -- so the "order by
        # similarity" toggle's client-side rebuild (see SHARED_JS's
        # SYN.buildDotplotLayout) reuses these rather than recomputing them
        'dpRx': dp_rx, 'dpRy': dp_ry, 'dpTotalX': dp_total_x, 'dpTotalY': dp_total_y,
        # None (-> JSON null, falsy in JS) if --stats wasn't given -- every
        # reader of this (SYN.buildStatsHtml) already treats that as "no
        # panel to show"
        'alignmentStats': alignment_stats,
    }
    # Seeds SYN.state's dotplot offsets to match Python's own initial render
    # (natural order -- see Dataset.__init__'s reversed dp_subject_offsets
    # and SHARED_JS's SYN.dpNaturalOrder) so the color/palette/min-block-
    # anchors callbacks have a valid SYN.state.dpQueryOffsets/dpSubjectOffsets
    # to redraw against before the order-by-similarity toggle is ever touched.
    # targetLabel/comparisonLabel seed to the same default the label
    # TextInputs show, so the very first pivot/pair click (before either
    # input is ever touched) already reads a real string, not null.
    seed_js = """
    (function() {
        const natural = SYN.dpNaturalOrder();
        SYN.state.dpQueryOffsets = SYN.computeOffsets(natural.queryOrder, SYN.data.querySizes);
        SYN.state.dpSubjectOffsets = SYN.computeOffsets(natural.subjectOrder, SYN.data.subjectSizes);
        SYN.state.targetLabel = SYN.data.targetName;
        SYN.state.comparisonLabel = SYN.data.comparisonName;
    })();
    """
    injected = ("<script>\n" + SHARED_JS + "\nSYN.data = " + json.dumps(syn_data) + ";\n"
                + seed_js + "\n</script>\n</body>")
    return page_html.replace("</body>", injected, 1), links_by_comparison, links_by_query


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--query_name', required=True)
    parser.add_argument('--subject_name', required=True)
    parser.add_argument('--query_chrom_sizes', required=True)
    parser.add_argument('--subject_chrom_sizes', required=True)
    parser.add_argument('--links', required=True)
    parser.add_argument('--target_homeolog_links', default=None,
                         help='optional target-vs-itself links TSV -- drawn as ribbons '
                              'within the target half of the ring (not included in the '
                              'zoom panel, which only shows target<->comparison links for the '
                              'clicked chromosome)')
    parser.add_argument('--comparison_homeolog_links', default=None,
                         help='optional comparison-vs-itself links TSV -- drawn as ribbons '
                              'within the comparison half of the ring (ditto, not included in '
                              'the zoom panel)')
    parser.add_argument('--stats', default=None,
                         help='optional compute_alignment_stats.py JSON -- shown as a small '
                              'always-visible summary panel on the page (proteome size, each '
                              "genome's aligned-protein count/mean identity) rather than baked "
                              'into any one exportable figure')
    parser.add_argument('--query_subtitle', default=None,
                         help='optional subtitle under the query/target name in the page '
                              'header -- e.g. the input file name, so a viewer can tell '
                              'which physical genome "target"/"comparison" refer to')
    parser.add_argument('--subject_subtitle', default=None,
                         help='ditto, under the subject/comparison name')
    parser.add_argument('--out_prefix', required=True)
    args = parser.parse_args()

    ds = Dataset(args.query_chrom_sizes, args.subject_chrom_sizes, args.links,
                 args.target_homeolog_links, args.comparison_homeolog_links)
    if not ds.query_chroms or not ds.subject_chroms:
        sys.exit("ERROR: no chromosomes to plot -- check --min_seq_size isn't "
                  "filtering out everything")

    alignment_stats = None
    if args.stats:
        with open(args.stats) as f:
            alignment_stats = json.load(f)

    page_html, links_by_comparison, links_by_query = build_page(
        ds, args.query_name, args.subject_name, args.query_subtitle, args.subject_subtitle,
        alignment_stats)

    out_path = f"{args.out_prefix}.interactive.html"
    with open(out_path, 'w') as f:
        f.write(page_html)

    total_links = sum(len(v) for v in links_by_comparison.values())
    homeolog_bits = []
    if args.target_homeolog_links:
        homeolog_bits.append(f"{len(ds.target_homeolog_links)} target homeolog link(s)")
    if args.comparison_homeolog_links:
        homeolog_bits.append(f"{len(ds.comparison_homeolog_links)} comparison homeolog link(s)")
    print(f"[plot_synteny_interactive] wrote {out_path} ({len(links_by_comparison)} comparison / "
          f"{len(links_by_query)} target chromosome(s), {total_links} link(s) embedded"
          + ((f"; {'; '.join(homeolog_bits)}") if homeolog_bits else '') + ")",
          file=sys.stderr)


if __name__ == '__main__':
    main()
