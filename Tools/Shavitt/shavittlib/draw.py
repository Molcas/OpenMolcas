#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2026, Roland Lindh                                     *
#***********************************************************************
"""Drawing of a Shavitt graph with matplotlib (no pyplot: usable from the GUI
and in batch mode alike).

Layout: the levels are the rows (head at the top, tail at the bottom); the
columns are the distinct (a,b) pairs, a descending, then b descending, so the
d=0 arcs are vertical and d=1,2,3 slant to the right. Highlighted walks are
drawn on top, each in its own colour slot of a fixed categorical order.
"""

import os
from dataclasses import dataclass, field

import numpy as np

# fixed categorical order (validated palette, light mode)
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
INK = '#1f1f1e'
INK2 = '#5f5e58'
GRID = '#c9c8c0'
FAINT = '#e4e3dc'
SURFACE = '#ffffff'
MAX_WALKS = len(SERIES)


@dataclass
class Walk:
    """A highlighted walk: label for the legend, step vector, weight c^2 (or None), colour slot 0..7."""
    label: str
    steps: np.ndarray
    weight: float = None
    slot: int = 0


BRA = '#6250d6'  # violet: bra (left) walk of a loop
KET = '#eb6834'  # orange, dashed: ket (right) walk of a loop


@dataclass
class LoopOverlay:
    """Loops to draw: the loops (shavittlib.loops.Loop), the arcs of the walk parts common to bra and
    ket (upper, lower), the selected midvertices, and whether each loop is drawn in detail."""
    loops: list = field(default_factory=list)
    context: list = field(default_factory=list)  # arcs (upper vertex, lower vertex)
    midvertices: list = field(default_factory=list)
    detail: bool = True
    note: str = ''


@dataclass
class SegmentOverlay:
    """Segments to draw (shavittlib.segments.Segment) and the vertex they belong to."""
    segments: list = field(default_factory=list)
    vertex: int = None
    detail: bool = True
    note: str = ''


@dataclass
class DrawOptions:
    vertex_numbers: bool = False  # number the vertices (1-based, as in the DRT)
    paldus_labels: bool = False  # (a,b,c) at each vertex
    midlevel_line: bool = True  # dashed line at the midlevel
    mark_midvertices: bool = False  # midvertices as squares
    step_labels: bool = False  # step d at the middle of each arc
    dim_graph: bool = True  # fade the graph behind highlighted walks
    width_by_weight: bool = False  # line width of a walk ~ sqrt(weight)
    level_orbitals: bool = False  # active orbital index in the level labels
    selected_vertex: int = None  # 0-based; ringed, with its subgraph if show_subgraph
    show_subgraph: bool = True
    show_logo: bool = True  # OpenMolcas logo in the lower right corner
    logo_color: str = 'black'  # any matplotlib color; the logo is drawn in it (one color, transparent background)
    title: str = ''
    extra: dict = field(default_factory=dict)


def layout(g):
    """x, y coordinates of all vertices and the columns in use as {x: (a,b)}.

    The columns form a regular (a,b) lattice, x = (amax-a)*(bmax+1) + (bmax-b),
    empty columns included, so that every step type has a fixed displacement
    (d=0: 0, d=1: 1, d=2: bmax, d=3: bmax+1) and repeated steps are straight lines.
    """
    a = g.drt[:, 2].astype(int)
    b = g.drt[:, 3].astype(int)
    bmax = int(b.max())
    x = (int(a.max()) - a) * (bmax + 1) + (bmax - b)
    x = (x - x.min()).astype(float)
    y = g.drt[:, 0].astype(float)
    cols = {int(xi): (int(ai), int(bi)) for xi, ai, bi in zip(x, a, b)}
    return x, y, cols


def figure_size(g):
    x, _, _ = layout(g)
    return max(6.0, 0.32 * (x.max() + 1) + 4.5), max(4.5, 0.62 * g.nlev + 2.2)


def draw(g, ax, walks=(), opts=None, loops=None, segments=None):
    """Draw the graph g on the matplotlib Axes ax. Returns the vertex coordinates (x, y).

    loops: a LoopOverlay of coupling-coefficient loops, or None.
    """
    if opts is None:
        opts = DrawOptions()
    x, y, cols = layout(g)
    ncol = int(x.max()) + 1
    walks = list(walks)[:MAX_WALKS]
    ax.clear()

    # midlevel
    if opts.midlevel_line:
        ax.axhline(g.midlev, color=GRID, lw=0.8, ls=(0, (4, 3)), zorder=0)
        ax.text(ncol - 0.5, g.midlev, ' midlevel', color=INK2, va='center', ha='left', fontsize=8)

    # all arcs
    base = FAINT if (walks and opts.dim_graph) else GRID
    for iv, ic, d in g.arcs():
        ax.plot([x[iv], x[ic]], [y[iv], y[ic]], color=base, lw=1.0, zorder=1, solid_capstyle='round')
        if opts.step_labels:
            xm, ym = (x[iv] + x[ic]) / 2, (y[iv] + y[ic]) / 2
            ax.text(xm, ym, str(d), fontsize=6, color=INK2, ha='center', va='center', zorder=2,
                    bbox=dict(boxstyle='round,pad=0.1', fc=SURFACE, ec='none'))

    # subgraph of the walks through the selected vertex
    sel = opts.selected_vertex
    if sel is not None and opts.show_subgraph:
        for iv, ic, d in g.arcs_through(sel):
            ax.plot([x[iv], x[ic]], [y[iv], y[ic]], color=INK2, lw=1.6, zorder=2, solid_capstyle='round')

    # highlighted walks, slightly offset so that shared arcs stay visible
    nw = len(walks)
    wmax = max([w.weight for w in walks if w.weight] or [1.0])
    for k, w in enumerate(walks):
        path = g.walk(w.steps)
        dx = (k - (nw - 1) / 2) * 0.06
        lw = 2.0
        if opts.width_by_weight and w.weight:
            lw = 1.0 + 4.0 * np.sqrt(w.weight / wmax)
        ax.plot(x[path] + dx, y[path], color=SERIES[w.slot % MAX_WALKS], lw=lw, zorder=3 + k,
                solid_capstyle='round', solid_joinstyle='round', label=w.label)

    if loops is not None:
        _draw_loops(g, ax, x, y, loops)
    if segments is not None:
        _draw_segments(g, ax, x, y, segments)

    # vertices with a surface ring; midvertices as squares
    mid = set(g.midvertices()) if opts.mark_midvertices else set()
    ax.scatter(x, y, s=46, color=INK, edgecolors=SURFACE, linewidths=1.5, zorder=20)
    if mid:
        m = sorted(mid)
        ax.scatter(x[m], y[m], s=150, marker='s', facecolors='none', edgecolors=INK, linewidths=1.4, zorder=21,
                   label='midvertex')
    if sel is not None:
        ax.scatter([x[sel]], [y[sel]], s=260, facecolors='none', edgecolors=INK, linewidths=1.6, zorder=21)

    for iv in range(g.nvert):
        if opts.vertex_numbers:
            ax.annotate(f'{iv + 1}', (x[iv], y[iv]), xytext=(5, 4), textcoords='offset points', fontsize=6.5,
                        color=INK, zorder=22)
        if opts.paldus_labels:
            a, b, c = g.drt[iv, 2:5]
            ax.annotate(f'{a},{b},{c}', (x[iv], y[iv]), xytext=(5, -4), textcoords='offset points', fontsize=6,
                        color=INK2, va='top', zorder=22)

    # axes: levels with their irreps on the left, (a,b) columns on top
    ax.set_yticks(range(g.nlev + 1))
    lab = ['0']
    for L in range(1, g.nlev + 1):
        t = f'{L}  {g.irrep_name(g.level_sym[L - 1])}'
        if opts.level_orbitals:
            t += f'  ({g.level_to_act[L - 1]})'
        lab.append(t)
    ax.set_yticklabels(lab)
    ax.set_ylabel('level  irrep' + ('  (orbital)' if opts.level_orbitals else ''), color=INK2)
    ax.set_xticks(sorted(cols))
    ax.set_xticklabels([f'{cols[k][0]},{cols[k][1]}' for k in sorted(cols)], fontsize=7, rotation=90)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.set_xlabel('(a,b)', color=INK2)
    ax.set_xlim(-0.6, ncol - 0.4 + 1.2)
    ax.set_ylim(-0.5, g.nlev + 0.5)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(colors=INK2, length=0)
    ax.set_facecolor(SURFACE)

    if opts.show_logo:
        add_logo(ax, color=opts.logo_color)

    if opts.title:
        ax.set_title(opts.title, color=INK, fontsize=10, loc='left', pad=30)
    if nw or mid or (loops is not None and loops.loops) or (segments is not None and segments.segments):
        ax.legend(title=WALKS_TITLE if nw else None, **LEGEND_STYLE)
    return x, y


LEGEND_STYLE = dict(loc='upper left', bbox_to_anchor=(1.0, 1.0), frameon=False, fontsize=8, labelcolor=INK,
                    title_fontsize=8)
WALKS_TITLE = 'CSF  step vector  coeff'
NOTE_LINES = 2  # the longest note under the graph (a loop value and a common walk)


def _widest_legend(g, kind):
    """The widest legend entries (and title) the overlay kind 'walks', 'loops' or 'segments' can show."""
    if kind == 'loops':
        return ['common walk parts', 'bra (L)', 'ket (R)', 'bra (L) arcs', 'ket (R) arcs', 'loop midvertex',
                'midvertex'], None
    if kind == 'segments':
        return ['bra (L) step', 'ket (R) step', 'segment vertex', 'midvertex'], None
    steps = '2' * g.nlev
    irr = max((g.irrep_name(i) for i in range(g.nsym)), key=len)
    return [f'{g.ncsf.max():>6d}  {steps}  -0.000', f'     -  {steps}  ({irr})', 'midvertex'], WALKS_TITLE


def fixed_layout(fig, ax, g, kind='walks'):
    """Lay out the figure with room for the widest legend and note of this overlay kind.

    Unlike tight_layout with the legend in it, the axes do not move when the legend or the
    note under the graph comes and goes (a vertex picked or unpicked).
    """
    from matplotlib.legend import Legend
    from matplotlib.lines import Line2D
    renderer = fig.canvas.get_renderer()
    labels, title = _widest_legend(g, kind)
    probe = Legend(ax, [Line2D([], [], lw=2.0) for _ in labels], labels, title=title, **LEGEND_STYLE)
    probe.set_figure(fig)
    w = probe.get_window_extent(renderer).width + 8
    h = 0.0
    if kind in ('loops', 'segments'):
        h = NOTE_LINES * 8 * 1.25 * fig.dpi / 72.0 + 6
    leg = ax.get_legend()
    if leg is not None:
        leg.set_in_layout(False)
    for t in ax.texts:
        if t.get_gid() == 'note':
            t.set_in_layout(False)
    fw, fh = fig.get_figwidth() * fig.dpi, fig.get_figheight() * fig.dpi
    fig.tight_layout(rect=(0.0, min(h / fh, 0.3), 1.0 - min(w / fw, 0.5), 1.0))


def _draw_loops(g, ax, x, y, ov):
    """The loop overlay: common walk parts dark, bra solid violet, ket dashed orange."""
    for u, c in ov.context:
        ax.plot([x[u], x[c]], [y[u], y[c]], color=INK2, lw=1.8, zorder=2.5, solid_capstyle='round')
    if ov.context:
        ax.plot([], [], color=INK2, lw=1.8, label='common walk parts')
    if ov.detail:
        for k, lp in enumerate(ov.loops):
            bx, by = x[lp.bra] - 0.07, y[lp.bra]
            kx, ky = x[lp.ket] + 0.07, y[lp.ket]
            ax.plot(bx, by, color=BRA, lw=2.6, zorder=12, solid_capstyle='round', solid_joinstyle='round',
                    label='bra (L)' if k == 0 else None)
            ax.plot(kx, ky, color=KET, lw=2.6, ls=(0, (3, 1.6)), zorder=13, dash_capstyle='round',
                    label='ket (R)' if k == 0 else None)
            if lp.j is not None and lp.kind in ('upper', 'lower', 'upper-partial'):
                v = lp.bra[0]
                ax.annotate(f'j={lp.j}', (x[v], y[v]), xytext=(-8, 6), textcoords='offset points', ha='right',
                            fontsize=8, color=INK, zorder=23, fontweight='bold')
            if lp.i is not None and lp.kind in ('upper', 'lower', 'lower-partial'):
                v = lp.bra[-1]
                ax.annotate(f'i={lp.i}', (x[v], y[v]), xytext=(-8, -12), textcoords='offset points', ha='right',
                            fontsize=8, color=INK, zorder=23, fontweight='bold')
    else:
        # envelope: every arc used by a bra or a ket walk of the loops, once
        ba, ka = set(), set()
        for lp in ov.loops:
            ba.update(zip(lp.bra[:-1], lp.bra[1:]))
            ka.update(zip(lp.ket[:-1], lp.ket[1:]))
        for k, (u, c) in enumerate(sorted(ba)):
            ax.plot([x[u] - 0.07, x[c] - 0.07], [y[u], y[c]], color=BRA, lw=1.8, alpha=0.75, zorder=12,
                    solid_capstyle='round', label='bra (L) arcs' if k == 0 else None)
        for k, (u, c) in enumerate(sorted(ka)):
            ax.plot([x[u] + 0.07, x[c] + 0.07], [y[u], y[c]], color=KET, lw=1.8, alpha=0.75, zorder=13,
                    ls=(0, (3, 1.6)), label='ket (R) arcs' if k == 0 else None)
    if ov.midvertices:
        m = ov.midvertices
        ax.scatter(x[m], y[m], s=150, marker='s', color=INK, edgecolors=SURFACE, linewidths=1.5, zorder=24,
                   label='loop midvertex')
    if ov.note:
        ax.text(0.0, -0.02, ov.note, transform=ax.transAxes, fontsize=8, color=INK, va='top', ha='left', gid='note')


def _draw_segments(g, ax, x, y, ov):
    """The segment overlay: bra arcs solid violet, ket arcs dashed orange, the vertex ringed."""
    lw, alpha = (2.8, 1.0) if ov.detail else (1.6, 0.6)
    seen_b, seen_k = set(), set()
    for s in ov.segments:
        same = s.vlt == s.vrt and s.vlb == s.vrb
        dx = 0.0 if same else 0.07
        if (s.vlt, s.vlb) not in seen_b:
            seen_b.add((s.vlt, s.vlb))
            ax.plot([x[s.vlt] - dx, x[s.vlb] - dx], [y[s.vlt], y[s.vlb]], color=BRA, lw=lw, alpha=alpha, zorder=12,
                    solid_capstyle='round', label='bra (L) step' if len(seen_b) == 1 else None)
        if (s.vrt, s.vrb) not in seen_k:
            seen_k.add((s.vrt, s.vrb))
            ax.plot([x[s.vrt] + dx, x[s.vrb] + dx], [y[s.vrt], y[s.vrb]], color=KET, lw=lw, alpha=alpha, zorder=13,
                    ls=(0, (3, 1.6)), label='ket (R) step' if len(seen_k) == 1 else None)
        if ov.detail:
            # the four corners of the segment
            pts = sorted({s.vlt, s.vrt, s.vlb, s.vrb})
            ax.scatter(x[pts], y[pts], s=70, facecolors='none', edgecolors=INK, linewidths=1.2, zorder=22)
            xm = (x[s.vlt] + x[s.vlb] + x[s.vrt] + x[s.vrb]) / 4
            ym = (y[s.vlt] + y[s.vlb]) / 2
            ax.annotate(f'#{s.number}  {s.value:+.4f}', (xm, ym), xytext=(10, 0), textcoords='offset points',
                        fontsize=8, color=INK, zorder=23, fontweight='bold', va='center',
                        bbox=dict(boxstyle='round,pad=0.2', fc=SURFACE, ec='none', alpha=0.85))
    if ov.vertex is not None:
        v = ov.vertex
        ax.scatter([x[v]], [y[v]], s=260, facecolors='none', edgecolors=INK, linewidths=1.8, zorder=24,
                   label='segment vertex')
    if ov.note:
        ax.text(0.0, -0.02, ov.note, transform=ax.transAxes, fontsize=8, color=INK, va='top', ha='left', gid='note')


LOGO = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'openmolcas_logo.png')


def add_logo(ax, height=0.11, margin=4, color='black'):
    """The OpenMolcas logo in the lower right corner of the figure.

    height: fraction of the figure height; margin: distance from the edges in points;
    color: any matplotlib color (the logo is one color on a transparent background).
    The logo belongs to ax (removed by ax.clear) but is placed in figure coordinates
    and kept out of the layout computation.
    """
    if not os.path.isfile(LOGO):
        return
    import matplotlib.image as mpimg
    from matplotlib.offsetbox import AnnotationBbox, OffsetImage
    from matplotlib.colors import to_rgb
    img = mpimg.imread(LOGO).copy()  # RGBA, values 0..1
    img[..., :3] = to_rgb(color)
    zoom = height * ax.figure.get_figheight() * 72.0 / img.shape[0]
    ab = AnnotationBbox(OffsetImage(img, zoom=zoom), (1.0, 0.0), xycoords='figure fraction',
                        xybox=(-margin, margin), boxcoords='offset points', box_alignment=(1.0, 0.0),
                        frameon=False, pad=0.0, annotation_clip=False, zorder=30)
    ab.set_in_layout(False)
    ax.add_artist(ab)


def save(g, filename, walks=(), opts=None, dpi=300, limits=None, loops=None, segments=None):
    """Draw the graph into a new figure and write it to filename (format from the extension).

    limits: ((xmin, xmax), (ymin, ymax)) to save a zoomed view.
    """
    from matplotlib.figure import Figure
    fig = Figure(figsize=figure_size(g), facecolor=SURFACE)
    ax = fig.add_subplot(111)
    draw(g, ax, walks, opts, loops, segments)
    if limits is not None:
        ax.set_xlim(limits[0])
        ax.set_ylim(limits[1])
    fig.tight_layout()
    # no bbox_inches='tight': it would crop the logo in the figure corner
    fig.savefig(filename, dpi=dpi, facecolor=SURFACE)
