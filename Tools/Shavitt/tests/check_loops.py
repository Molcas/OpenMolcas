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
"""Check the loops of shavittlib.loops.

usage: check_loops.py rasscf.h5 [--no-pairs]

1. The active one-body density matrix built from the loop values and the CI
   vector equals DENSITY_MATRIX of the file (the coupling coefficients).
2. The split-graph loops enumerated at the midvertices (upper and lower loops
   at one midvertex, upper and lower partial loops at a pair) are exactly the
   loops found between all pairs of walks of the graph, with the same values
   (skipped with --no-pairs; the cost grows with the square of the walks).
3. The common walks of every loop (CommonWalks: one complete loop, or one
   upper and one lower partial loop) are distinct pairs of complete walks
   whose coupling coefficient is the loop value, and together they cover the
   common parts drawn by the GUI. With pairs: all the embeddings together are
   exactly the pairs of walks connected by an E_ij.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from shavittlib.graph import ShavittGraph  # noqa: E402
from shavittlib.loops import (CommonWalks, _walk_pair_value, common_parts, coupling,  # noqa: E402
                              crossing_loops, density_matrix, midvertex_loops, midvertex_partners)


def all_walks(g):
    out = []

    def rec(v, L, steps):
        if L == 0:
            out.append(list(reversed(steps)))
            return
        for d in range(4):
            c = g.down[v, d]
            if c >= 0:
                rec(int(c), L - 1, steps + [d])
    rec(0, g.nlev, [])
    return out


def main(fn, pairs=True):
    import h5py
    g = ShavittGraph.from_h5(fn)
    nbad = 0

    with h5py.File(fn, 'r') as f:
        Dr = np.asarray(f['DENSITY_MATRIX'][()])
        a2l = np.asarray(f.attrs['A2LEV'][:g.nlev]) - 1
    for root in range(g.ci_vectors.shape[0]):
        D = density_matrix(g, root)[np.ix_(a2l, a2l)]
        err = np.abs(D - Dr[root]).max()
        print(f'root {root + 1}: max |D(loops) - DENSITY_MATRIX| = {err:.2e}')
        if err > 1e-8:
            nbad += 1

    # 3. common walks: each embedding of a loop is a pair of walks with the loop value
    embedded = set()

    def check_common(sel, mv, i, j, value):
        bad = 0
        cw = CommonWalks(g, sel, mv)
        union = set()
        for lp in sel:
            union |= set(common_parts(g, lp, mv))
        seen = set()
        for k in range(cw.count):
            seen |= {(u, c) for u, c, _ in cw.arcs(k)}
            bra, ket = cw.step_vectors(k)
            key = (tuple(bra), tuple(ket))
            r = coupling(g, bra, ket)
            if key in embedded or r is None or r[3] or (r[0], r[1]) != (i, j) or abs(r[2] - value) > 1e-10:
                bad += 1
            embedded.add(key)
        return bad + (seen != union)

    nloop = 0
    for mv in g.midvertices():
        for lp in midvertex_loops(g, mv, limit=10**9):
            nbad += check_common([lp], mv, lp.i, lp.j, lp.value)
            nloop += 1
        for w in midvertex_partners(g, mv):
            if w < mv:
                continue
            _, _, ups, lows = crossing_loops(g, mv, w, limit=10**9)
            for u in ups:
                for lo in lows:
                    nbad += check_common([u, lo], None, lo.i, u.j, u.value * lo.value)
                    nloop += 1
    print(f'common walks: {nloop} loops embedded in {len(embedded)} pairs of walks')

    if pairs:
        m = g.midlev
        walks = all_walks(g)
        paths = [g.walk(np.array(s)) for s in walks]

        def vx(p, L):
            return p[g.nlev - L]
        expect = set()
        for a in range(len(walks)):
            for b in range(len(walks)):
                if a == b:
                    continue
                r = _walk_pair_value(g, walks[a], walks[b])
                if r is None:
                    continue
                i, j, v = r
                pl, pr = paths[a], paths[b]
                if i > m or j <= m:
                    kind = 'upper' if i > m else 'lower'
                    mv = vx(pl, m)
                    expect.add((kind, mv, tuple(vx(pl, L) for L in range(j, i - 2, -1)),
                                tuple(vx(pr, L) for L in range(j, i - 2, -1)), round(v, 10)))
                else:
                    expect.add(('upper-partial', None, tuple(vx(pl, L) for L in range(j, m - 1, -1)),
                                tuple(vx(pr, L) for L in range(j, m - 1, -1)), None))
                    expect.add(('lower-partial', None, tuple(vx(pl, L) for L in range(m, i - 2, -1)),
                                tuple(vx(pr, L) for L in range(m, i - 2, -1)), None))
                    expect.add(('crossing', None, (vx(pl, m), vx(pr, m)), (j, i), None))
        found = set()
        prods = {}
        for mv in g.midvertices():
            for lp in midvertex_loops(g, mv, limit=10**9):
                found.add((lp.kind, mv, tuple(lp.bra), tuple(lp.ket), round(lp.value, 10)))
            for w in midvertex_partners(g, mv):
                L, R, ups, lows = crossing_loops(g, mv, w, limit=10**9)
                for u in ups:
                    found.add(('upper-partial', None, tuple(u.bra), tuple(u.ket), None))
                for lo in lows:
                    found.add(('lower-partial', None, tuple(lo.bra), tuple(lo.ket), None))
                for u in ups:
                    for lo in lows:
                        found.add(('crossing', None, (L, R), (u.j, lo.i), None))
                        prods[(tuple(u.bra[:-1] + lo.bra), tuple(u.ket[:-1] + lo.ket))] = u.value * lo.value
        # values of the crossing loops: product of the partial loops
        for a in range(len(walks)):
            for b in range(len(walks)):
                r = _walk_pair_value(g, walks[a], walks[b]) if a != b else None
                if r is None or not (r[0] <= m < r[1]):
                    continue
                i, j, v = r
                key = (tuple(paths[a][g.nlev - L] for L in range(j, i - 2, -1)),
                       tuple(paths[b][g.nlev - L] for L in range(j, i - 2, -1)))
                if key not in prods or abs(prods[key] - v) > 1e-10:
                    nbad += 1
                    print('crossing loop value mismatch', key, v, prods.get(key))
        connected = {(tuple(walks[a]), tuple(walks[b])) for a in range(len(walks)) for b in range(len(walks))
                     if a != b and _walk_pair_value(g, walks[a], walks[b]) is not None}
        if connected != embedded:
            nbad += 1
            print(f'embeddings {len(embedded)} differ from the connected pairs of walks {len(connected)}')
        missing, extra = expect - found, found - expect
        print(f'{len(walks)} walks; loops by pairs of walks {len(expect)}, enumerated {len(found)}; '
              f'missing {len(missing)}, extra {len(extra)}')
        for x in list(missing)[:5]:
            print('  missing', x)
        for x in list(extra)[:5]:
            print('  extra', x)
        nbad += len(missing) + len(extra)
    print('OK' if nbad == 0 else 'FAILED')
    return 1 if nbad else 0


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1], '--no-pairs' not in sys.argv))
