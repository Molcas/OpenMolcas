#!/usr/bin/env python3
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
"""Shavitt: the Shavitt graph of a RASSCF/CASSCF wave function.

Without -o or --text a graphical interface opens; the display options given
on the command line are its initial settings.

examples:
  shavitt.py h2o.rasscf.h5                         # graphical interface
  shavitt.py                                       # ... and open a file from it
  shavitt.py h2o.rasscf.h5 --text                  # summary and leading CSFs
  shavitt.py h2o.rasscf.h5 -o h2o.pdf              # draw root 1 to a file
  shavitt.py h2o.rasscf.h5 -r 2 -t 0.1 -n 8 -o g.png --numbers --midvertices
  shavitt.py h2o.rasscf.h5 --csf 18 27 -o g.pdf    # highlight given CSFs
  shavitt.py h2o.rasscf.h5 --walk 2ud2ud0 -o g.pdf # highlight a step vector

a graph from scratch (no file; orbitals per irrep, D2h subgroups):
  shavitt.py --electrons 8 --cas 3 1 2 1           # CAS(8,7), 4 irreps, singlet A1
  shavitt.py --electrons 8 --spin 3 --irrep 2 --ras1 1 0 1 0 --ras2 1 1 0 1 \
             --ras3 2 0 1 1 --holes 2 --particles 2 --text
  shavitt.py --electrons 6 --gas 2:2:4 --gas 2:4:6 --gas 2:6:6 -o gas.pdf
             # GAS: ORBITALS:MIN:MAX per space, the electrons of the spaces 1..k together
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shavittlib import __version__  # noqa: E402
from shavittlib.draw import MAX_WALKS, DrawOptions, Walk  # noqa: E402
from shavittlib.graph import ShavittGraph  # noqa: E402


def options(a):
    return DrawOptions(vertex_numbers=a.numbers, paldus_labels=a.paldus, midlevel_line=not a.no_midlevel,
                       mark_midvertices=a.midvertices, step_labels=a.steps, level_orbitals=a.orbitals,
                       width_by_weight=a.width_by_weight, dim_graph=not a.no_dim, show_logo=not a.no_logo, logo_color=a.logo_color,
                       selected_vertex=None if a.vertex is None else a.vertex - 1)


def graph_from_args(a):
    """The graph defined on the command line (--cas, --ras1..3 or --gas)."""
    if a.electrons is None:
        raise ValueError('--electrons is needed for a graph from scratch')
    gas = None
    if a.gas:
        gas = []
        for t in a.gas:
            try:
                orbs, lo, hi = t.split(':')
                gas.append(([int(x) for x in orbs.replace(',', ' ').split()], int(lo), int(hi)))
            except ValueError:
                raise ValueError(f'--gas {t}: expected ORBITALS:MIN:MAX, e.g. 2,0,1,0:2:4') from None
    return ShavittGraph.from_definition(a.electrons, a.spin, a.irrep, cas=a.cas, ras1=a.ras1, ras2=a.ras2, ras3=a.ras3,
                                        holes=a.holes, particles=a.particles, gas=gas, irrep_labels=a.irrep_labels)


def main(argv=None):
    p = argparse.ArgumentParser(prog='shavitt', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('h5', nargs='?', help='rasscf.h5 file')
    p.add_argument('--version', action='version', version=f'Shavitt {__version__}')
    p.add_argument('-r', '--root', type=int, default=1, help='root (1-based) whose CSFs are shown (default 1)')
    p.add_argument('-t', '--threshold', type=float, default=0.05, help='smallest |coefficient| shown (default 0.05)')
    p.add_argument('-n', '--max-walks', type=int, default=5, help=f'leading CSFs drawn (default 5, at most {MAX_WALKS})')
    p.add_argument('--csf', type=int, nargs='+', metavar='I', help='highlight these CSFs (1-based) instead')
    p.add_argument('--walk', nargs='+', metavar='STEPS', help='highlight these step vectors (0/u/d/2, level 1 first)')
    p.add_argument('-o', '--output', help='draw the graph to this file (pdf, png, svg, ...) without the interface')
    p.add_argument('--dpi', type=int, default=300, help='resolution of bitmap output (default 300)')
    p.add_argument('--text', action='store_true', help='print the summary and the CSFs only')
    p.add_argument('--vertices', action='store_true', help='print the distinct-row table')
    b_ = p.add_argument_group('graph from scratch (instead of a file)')
    b_.add_argument('--electrons', type=int, metavar='N', help='number of active electrons')
    b_.add_argument('--spin', type=int, default=1, metavar='2S+1', help='spin multiplicity (default 1)')
    b_.add_argument('--irrep', type=int, default=1, metavar='I', help='target irrep, 1-based (default 1)')
    b_.add_argument('--cas', type=int, nargs='+', metavar='N', help='CAS: active orbitals of each irrep')
    for k in (1, 2, 3):
        b_.add_argument(f'--ras{k}', type=int, nargs='+', metavar='N', help=f'RAS{k} orbitals of each irrep')
    b_.add_argument('--holes', type=int, default=0, help='RAS: maximum holes in RAS1')
    b_.add_argument('--particles', type=int, default=0, help='RAS: maximum electrons in RAS3')
    b_.add_argument('--gas', action='append', metavar='ORBS:MIN:MAX',
                    help='GAS space: orbitals of each irrep (comma separated), minimum and maximum number of '
                         'electrons in the spaces 1..k together; repeat for each space')
    b_.add_argument('--irrep-labels', nargs='+', metavar='LABEL', help='names of the irreps')
    g_ = p.add_argument_group('display')
    g_.add_argument('--numbers', action='store_true', help='number the vertices')
    g_.add_argument('--paldus', action='store_true', help='(a,b,c) at each vertex')
    g_.add_argument('--midvertices', action='store_true', help='mark the midvertices')
    g_.add_argument('--no-midlevel', action='store_true', help='no midlevel line')
    g_.add_argument('--steps', action='store_true', help='label the arcs with their step d')
    g_.add_argument('--orbitals', action='store_true', help='active orbital in the level labels')
    g_.add_argument('--width-by-weight', action='store_true', help='line width of a walk ~ |c|')
    g_.add_argument('--no-dim', action='store_true', help='do not fade the graph behind the walks')
    g_.add_argument('--no-logo', action='store_true', help='leave out the OpenMolcas logo')
    g_.add_argument('--logo-color', default='black', metavar='COLOR',
                    help='color of the OpenMolcas logo (a name or #rrggbb; default black)')
    g_.add_argument('--loops', type=int, nargs='+', metavar='V',
                    help='the coupling-coefficient loops at midvertex V (1-based vertex number), or at the '
                         'partner pair V W (crossing loops); with --loop-head J / --loop-tail I to filter')
    g_.add_argument('--loop-head', type=int, metavar='J', help='only loops with head level J')
    g_.add_argument('--loop-tail', type=int, metavar='I', help='only loops with tail level I')
    g_.add_argument('--vertex', type=int, metavar='V', help='select vertex V (1-based) and show its subgraph')
    a = p.parse_args(argv)

    g = None
    if a.cas or a.ras1 or a.ras2 or a.ras3 or a.gas:
        if a.h5:
            p.error('give either a rasscf.h5 file or the definition of a graph, not both')
        try:
            g = graph_from_args(a)
        except ValueError as e:
            p.error(str(e))

    if not (a.output or a.text or a.vertices):
        from shavittlib.gui import run
        return run(a.h5, a.root, a.threshold, options(a), graph=g)

    if g is None:
        if not a.h5:
            p.error('a rasscf.h5 file (or the definition of a graph) is needed with -o, --text or --vertices')
        try:
            g = ShavittGraph.from_h5(a.h5)
        except (OSError, KeyError, ValueError) as e:
            sys.exit(f'shavitt.py: {e}')
        for w in g.warnings:
            print(f'Warning: {w}', file=sys.stderr)
    print(g.summary())

    if a.vertices:
        print('\n vertex   L   N   a   b   c    down(d=0..3)')
        for v in g.vertices():
            dn = ' '.join(f'{c + 1:3d}' if c >= 0 else '  .' for c in v.down)
            print(f'{v.index + 1:7d} {v.level:3d} {v.nel:3d} {v.a:3d} {v.b:3d} {v.c:3d}   {dn}')

    root = a.root - 1
    rows = []
    if a.csf:
        for i in a.csf:
            c = None if g.ci_vectors is None else float(g.ci_vectors[root, i - 1])
            rows.append((i - 1, c, g.step_vector(i - 1)))
    elif g.ci_vectors is not None and not a.walk:
        rows = g.leading_csfs(root, a.threshold)[:max(1, a.max_walks) if a.output else None]

    if rows:
        e = ''
        if g.root_energies is not None and a.root <= len(g.root_energies):
            e = f', E = {g.root_energies[root]:.8f}'
        print(f'\nroot {a.root}{e}')
        print('     CSF  step vector (levels 1..nlev)       coeff   weight')
        for i, c, s in rows:
            cs = '' if c is None else f'{c:10.5f} {c * c:8.5f}'
            print(f'{i + 1:8d}  {g.csf_string(s):30s} {cs}')

    if a.output:
        from shavittlib.draw import save
        # with --loops the leading CSFs are not drawn (as on the Loops tab), unless asked for with --csf
        drawn = rows if (a.csf or not a.loops) else []
        walks = [Walk(f'{i + 1:>6d}  {g.csf_string(s)}' + ('' if c is None else f'  {c:+.3f}'), s,
                      None if c is None else c * c, k) for k, (i, c, s) in enumerate(drawn[:MAX_WALKS])]
        for t in a.walk or []:
            s = g.parse_steps(t)
            i = g.csf_index(s)
            lab = (f'{i + 1:>6d}  {g.csf_string(s)}' if i is not None
                   else f'     -  {g.csf_string(s)}  ({g.irrep_name(g.walk_irrep(s))})')
            if len(walks) < MAX_WALKS:
                walks.append(Walk(lab, s, None, len(walks)))
        opts = options(a)
        opts.title = f'{g.ci_type} {g.nactel}e/{g.nlev}o, 2S+1={g.ispin}, irrep {g.irrep_name(g.lsym)}, root {a.root}'
        overlay = None
        if a.loops:
            from shavittlib import loops as lps
            from shavittlib.draw import LoopOverlay
            mvs = [v - 1 for v in a.loops[:2]]
            if any(v not in g.midvertices() for v in mvs):
                p.error(f'--loops: the midvertices are the vertices {[v + 1 for v in g.midvertices()]}')
            if len(mvs) == 1:
                items = lps.midvertex_loops(g, mvs[0], a.loop_head, a.loop_tail)
                print(f'\n{len(items)} loops at midvertex {a.loops[0]}')
            else:
                L, R, ups, lows = lps.crossing_loops(g, mvs[0], mvs[1], a.loop_head, a.loop_tail)
                items = ups + lows
                print(f'\ncrossing loops at L = vertex {L + 1}, R = vertex {R + 1}: '
                      f'{len(ups)} upper x {len(lows)} lower partial loops')
            for lp in items:
                print(f'  {lp.kind:13s} {lp.operator():9s} {lp.value:+.5f}  L {lp.steps_string("bra")}  '
                      f'R {lp.steps_string("ket")}')
            overlay = LoopOverlay(loops=items, midvertices=mvs, detail=False)
        save(g, a.output, walks, opts, dpi=a.dpi, loops=overlay)
        print(f'\ngraph written to {a.output}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
