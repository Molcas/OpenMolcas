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
"""Loops of the one-body coupling coefficients <L|E_ij|R> on a Shavitt graph,
and their split-graph partition at the midlevel.

The segment catalog and values are those of SGUGA (src/sguga_util/sguga.F90,
MKSEG; the ASTA segment values), segments 1..26 of the singlet operator E_ij:

  1- 4  walk segments above the loop (bra and ket on the same vertex)
  5- 8  head segments (at the upper level j of the loop)
  9-18  intermediate segments, db = -1 (9-13) and db = +1 (14-18)
 19-22  tail segments (at the lower level i of the loop)
 23-26  walk segments below the loop

The bra (left) walk L and the ket (right) walk R differ between the levels
i < j. Inside the loop the vertex pair (L, R) is of class 1 (a_R = a_L,
b_R = b_L - 1) or 2 (a_R = a_L - 1, b_R = b_L + 1): R has one electron less,
so L = E_ij R up to the coupling coefficient (an electron of level j in R is
at level i in L). The value of a loop is the product of its segment values.

Split graph: with the midlevel m, a loop lies in the upper half (i > m), in
the lower half (j <= m), or it crosses the midlevel (i <= m < j). A crossing
loop passes the midlevel at a pair of midvertices (L, R) of class 1 or 2 and
is the product of an upper and a lower partial loop.
"""

from dataclasses import dataclass, field

import numpy as np

# segments 1..26 of SGUGA (index 0..25 here)
ITVPT = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3]
IBVPT = [0, 0, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 2, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3]
IC1 = [0, 1, 2, 3, 0, 2, 0, 1, 0, 1, 1, 2, 3, 0, 1, 2, 2, 3, 1, 3, 2, 3, 0, 1, 2, 3]
IC2 = [0, 1, 2, 3, 1, 3, 2, 3, 0, 1, 2, 2, 3, 0, 1, 1, 2, 3, 0, 2, 0, 1, 0, 1, 2, 3]
ISVC = [1, 1, 1, 1, 1, 7, 8, 4, 1, 2, 9, 10, 2, 1, 2, 11, 12, 2, 1, 5, 6, 3, 1, 1, 1, 1]
HEAD = range(4, 8)
MID = range(8, 18)
TAIL = range(18, 22)


def segment_value(isvc, b):
    """The ASTA segment value of formula isvc for b of the top left vertex (MKSEG)."""
    s = -1.0 if b % 2 else 1.0
    if isvc == 0:
        return 0.0
    if isvc == 1:
        return 1.0
    if isvc == 2:
        return -1.0
    if isvc == 3:
        return s
    if isvc == 4:
        return -s
    if isvc == 5:
        return np.sqrt(b / (1.0 + b))
    if isvc == 6:
        return s * np.sqrt((1.0 + b) / (2.0 + b))
    if isvc == 7:
        return np.sqrt((2.0 + b) / (1.0 + b))
    if isvc == 8:
        return s * np.sqrt((2.0 + b) / (1.0 + b))
    if isvc == 9:
        return s / np.sqrt(b * (1.0 + b)) if b > 0 else 0.0
    if isvc == 10:
        return np.sqrt(b * (2.0 + b)) / (1.0 + b)
    if isvc == 11:
        return s / np.sqrt((1.0 + b) * (2.0 + b))
    if isvc == 12:
        return np.sqrt((1.0 + b) * (3.0 + b)) / (2.0 + b)
    if isvc == 13:
        return 2.0
    raise ValueError(isvc)


# (top class, bottom class, bra step, ket step) -> segment index, for the loop segments 5..22
_SEG = {}
for _k in list(HEAD) + list(MID) + list(TAIL):
    _SEG[(ITVPT[_k], IBVPT[_k], IC1[_k], IC2[_k])] = _k


def pair_class(g, vl, vr):
    """Class of the vertex pair (L, R) at one level: 0 equal, 1/2 loop interior, None otherwise."""
    if vl == vr:
        return 0
    al, bl = g.drt[vl, 2], g.drt[vl, 3]
    ar, br = g.drt[vr, 2], g.drt[vr, 3]
    if g.drt[vl, 0] != g.drt[vr, 0]:
        return None
    if ar == al and br == bl - 1:
        return 1
    if ar == al - 1 and br == bl + 1:
        return 2
    return None


@dataclass
class Loop:
    """A (partial) loop: the bra/ket vertices from level top down to level bottom-1 and their steps.

    kind: 'upper' or 'lower' (complete loop in one half), 'upper-partial' or
    'lower-partial' (the two halves of a loop crossing the midlevel), 'full'.
    j, i: head and tail level of the loop (None for the end that lies in the other half).
    """
    kind: str
    top: int  # level of the first segment
    bottom: int  # level of the last segment
    bra: list  # vertices, levels top .. bottom-1
    ket: list
    bra_steps: list  # steps, levels top .. bottom
    ket_steps: list
    value: float
    j: int = None
    i: int = None
    extra: dict = field(default_factory=dict)

    def steps_string(self, which):
        s = self.bra_steps if which == 'bra' else self.ket_steps
        return ''.join('0ud2'[d] for d in reversed(s))  # lowest level first, as a step vector

    def operator(self):
        a = '?' if self.i is None else str(self.i)
        b = '?' if self.j is None else str(self.j)
        return f'E({a},{b})'


def _descend(g, vl, vr, cls, level, stop, allow_tail, allow_mid_stop, out, segs, val, bra, ket):
    """Depth-first extension of a loop below level (inclusive) with interior class cls."""
    for k in list(MID) + (list(TAIL) if allow_tail(level) else []):
        if ITVPT[k] != cls:
            continue
        cl, cr = IC1[k], IC2[k]
        bl, br = g.down[vl, cl], g.down[vr, cr]
        if bl < 0 or br < 0:
            continue
        cb = pair_class(g, bl, br)
        if k in TAIL:
            if bl != br:
                continue
        elif cb != IBVPT[k]:
            continue
        v = val * segment_value(ISVC[k], int(g.drt[vl, 3]))
        if v == 0.0:
            continue
        nb, nk, ns = bra + [bl], ket + [br], segs + [(cl, cr)]
        if k in TAIL:
            out.append((level, nb, nk, ns, v, True))
        elif level - 1 == stop:
            if allow_mid_stop:
                out.append((level, nb, nk, ns, v, False))
        else:
            _descend(g, bl, br, cb, level - 1, stop, allow_tail, allow_mid_stop, out, ns, v, nb, nk)


def _loops_from_head(g, v, j, stop, allow_tail, allow_mid_stop):
    """All loops whose head is the vertex v at level j, ending at a tail above level stop, or at stop."""
    out = []
    for k in HEAD:
        cl, cr = IC1[k], IC2[k]
        bl, br = g.down[v, cl], g.down[v, cr]
        if bl < 0 or br < 0:
            continue
        cb = pair_class(g, bl, br)
        if cb != IBVPT[k]:
            continue
        val = segment_value(ISVC[k], int(g.drt[v, 3]))
        if val == 0.0:
            continue
        if j - 1 == stop:
            if allow_mid_stop:
                out.append((j, [v, bl], [v, br], [(cl, cr)], val, False))
            continue
        _descend(g, bl, br, cb, j - 1, stop, allow_tail, allow_mid_stop, out, [(cl, cr)], val, [v, bl], [v, br])
    return out


def _reach_down(g, v):
    """Vertices reachable downwards from v (v included)."""
    seen = {v}
    stack = [v]
    while stack:
        u = stack.pop()
        for c in g.down[u]:
            if c >= 0 and c not in seen:
                seen.add(int(c))
                stack.append(int(c))
    return seen


def _reach_up(g, v):
    """Vertices from which v is reachable (v included)."""
    up = {v}
    for L in range(int(g.drt[v, 0]) + 1, g.nlev + 1):
        for u in g.vertices_at_level(L):
            if any(c in up for c in g.down[u] if c >= 0):
                up.add(u)
    return up


def _mk(kind, j, rows):
    loops = []
    for (last, bra, ket, segs, v, closed) in rows:
        lp = Loop(kind, j, last, bra, ket, [s[0] for s in segs], [s[1] for s in segs], v)
        loops.append(lp)
    return loops


def midvertex_loops(g, mv, head=None, tail=None, upper=True, lower=True, limit=5000):
    """The loops of the split graph at one midvertex mv (0-based vertex index).

    upper: loops between the graph head and mv (both walks pass through mv);
    lower: loops between mv and the graph tail. Optional filters: head level j,
    tail level i. At most limit loops are returned.
    """
    m = g.midlev
    res = []
    if upper:
        above = _reach_up(g, mv)
        for j in range(g.nlev, m + 1, -1):
            if head is not None and j != head:
                continue
            for v in g.vertices_at_level(j):
                if v not in above:
                    continue
                rows = _loops_from_head(g, v, j, m, lambda L: L > m, False)
                for r in rows:
                    i = r[0]
                    if r[1][-1] not in above or (tail is not None and i != tail):
                        continue
                    lp = _mk('upper', j, [r])[0]
                    lp.j, lp.i = j, i
                    res.append(lp)
                    if len(res) >= limit:
                        return res
    if lower:
        below = _reach_down(g, mv)
        for j in range(m, 1, -1):
            if head is not None and j != head:
                continue
            for v in g.vertices_at_level(j):
                if v not in below:
                    continue
                rows = _loops_from_head(g, v, j, 0, lambda L: True, False)
                for r in rows:
                    i = r[0]
                    if tail is not None and i != tail:
                        continue
                    lp = _mk('lower', j, [r])[0]
                    lp.j, lp.i = j, i
                    res.append(lp)
                    if len(res) >= limit:
                        return res
    return res


def partner_class(g, mv1, mv2):
    """(L, R, class) for two midvertices that can carry a loop crossing the midlevel, else None."""
    c = pair_class(g, mv1, mv2)
    if c in (1, 2):
        return mv1, mv2, c
    c = pair_class(g, mv2, mv1)
    if c in (1, 2):
        return mv2, mv1, c
    return None


def midvertex_partners(g, mv):
    """The midvertices that form a crossing-loop pair with mv."""
    return [w for w in g.midvertices() if w != mv and partner_class(g, mv, w) is not None]


def crossing_loops(g, mv1, mv2, head=None, tail=None, limit=5000):
    """The upper and lower partial loops of the loops crossing the midlevel at the pair (mv1, mv2).

    Returns (L, R, upper partial loops, lower partial loops); every combination of an upper and a
    lower partial loop is a loop of <L..|E_ij|R..> with the value the product of the two.
    """
    pc = partner_class(g, mv1, mv2)
    if pc is None:
        raise ValueError('the midvertices are not a crossing-loop pair (same level, N differs by 1, b by 1)')
    L, R, c = pc
    m = g.midlev
    ups = []
    aboveL, aboveR = _reach_up(g, L), _reach_up(g, R)
    for j in range(g.nlev, m, -1):
        if head is not None and j != head:
            continue
        for v in g.vertices_at_level(j):
            if v not in aboveL or v not in aboveR:
                continue
            for r in _loops_from_head(g, v, j, m, lambda lev: False, True):
                if r[1][-1] == L and r[2][-1] == R:
                    lp = _mk('upper-partial', j, [r])[0]
                    lp.j = j
                    ups.append(lp)
                    if len(ups) >= limit:
                        break
    lows = []
    rows = []
    _descend(g, L, R, c, m, -1, lambda lev: True, False, rows, [], 1.0, [L], [R])
    for r in rows:
        i = r[0]
        if tail is not None and i != tail:
            continue
        lp = _mk('lower-partial', m, [r])[0]
        lp.i = i
        lows.append(lp)
        if len(lows) >= limit:
            break
    return L, R, ups, lows


def coupling(g, bra_steps, ket_steps):
    """The one-body coupling coefficient between two complete walks.

    Returns (i, j, value, swapped): <bra|E_ij|ket> (swapped False) or <ket|E_ij|bra> (swapped True),
    i < j levels; None if the walks are equal or not connected by a single E_ij.
    """
    for swapped, (sl, sr) in enumerate([(bra_steps, ket_steps), (ket_steps, bra_steps)]):
        r = _walk_pair_value(g, sl, sr)
        if r is not None:
            return r[0], r[1], r[2], bool(swapped)
    return None


def _walk_pair_value(g, sl, sr):
    pl, pr = g.walk(sl), g.walk(sr)  # vertices from the head: index k is level nlev-k
    diff = [L for L in range(1, g.nlev + 1) if pl[g.nlev - L] != pr[g.nlev - L] or sl[L - 1] != sr[L - 1]]
    if not diff:
        return None
    lv = [L for L in range(1, g.nlev + 1) if sl[L - 1] != sr[L - 1]]
    j, i = max(lv), min(lv)
    # vertices equal at and above level j, and at and below level i-1
    for L in range(j, g.nlev + 1):
        if pl[g.nlev - L] != pr[g.nlev - L]:
            return None
    for L in range(0, i):
        if pl[g.nlev - L] != pr[g.nlev - L]:
            return None
    val = 1.0
    for L in range(j, i - 1, -1):
        vl, vr = pl[g.nlev - L], pr[g.nlev - L]
        bl, br = pl[g.nlev - L + 1], pr[g.nlev - L + 1]
        ct = 0 if L == j else pair_class(g, vl, vr)
        cb = 3 if L == i else pair_class(g, bl, br)
        if ct is None or cb is None or (L == i and bl != br):
            return None
        k = _SEG.get((ct, cb, sl[L - 1], sr[L - 1]))
        if k is None:
            return None
        val *= segment_value(ISVC[k], int(g.drt[vl, 3]))
        if val == 0.0:
            return None
    return i, j, val


def density_matrix(g, root=0):
    """Active one-body density matrix (levels x levels) from the loops, brute force over CSF pairs (tests)."""
    n = int(g.ncsf[g.lsym])
    c = g.ci_vectors[root]
    steps = [g.step_vector(k) for k in range(n)]
    D = np.zeros((g.nlev, g.nlev))
    for k in range(n):
        for L in range(g.nlev):
            D[L, L] += c[k] ** 2 * (0, 1, 1, 2)[steps[k][L]]
    for a in range(n):
        for b in range(a + 1, n):
            r = _walk_pair_value(g, steps[a], steps[b])
            if r is None:
                r2 = _walk_pair_value(g, steps[b], steps[a])
                if r2 is None:
                    continue
                i, j, v = r2
            else:
                i, j, v = r
            D[i - 1, j - 1] += c[a] * c[b] * v
            D[j - 1, i - 1] += c[a] * c[b] * v
    return D


def arcs_between(g, top, bottom):
    """The arcs (upper, lower) of all walk pieces from vertex top down to vertex bottom."""
    if top == bottom:
        return []
    dn, up = _reach_down(g, top), _reach_up(g, bottom)
    return [(u, c) for u, c, d in g.arcs() if u in dn and u in up and c in up and c in dn]


def common_segments(g, lp, mv=None, upper=True, lower=True):
    """The pieces (top vertex, bottom vertex) of the walk that bra and ket of the loop lp share.

    upper/lower False: leave out the spectator half walk above/below the midvertex mv.
    """
    head, tail = 0, g.nvert - 1
    if lp.kind == 'upper':
        segs = [(head, lp.bra[0]), (lp.bra[-1], mv)] + ([(mv, tail)] if lower else [])
    elif lp.kind == 'lower':
        segs = ([(head, mv)] if upper else []) + [(mv, lp.bra[0]), (lp.bra[-1], tail)]
    elif lp.kind == 'upper-partial':
        segs = [(head, lp.bra[0])]
    elif lp.kind == 'lower-partial':
        segs = [(lp.bra[-1], tail)]
    else:
        segs = []
    return [(int(t), int(b)) for t, b in segs if t != b]


def common_parts(g, lp, mv=None, upper=True, lower=True):
    """The arcs of all the walk parts that bra and ket of the loop lp share (above and below the loop).

    upper/lower False: leave out the spectator half walk above/below the midvertex mv.
    """
    return [a for t, b in common_segments(g, lp, mv, upper, lower) for a in arcs_between(g, t, b)]


def walk_counts_to(g, bottom):
    """For every vertex v the number of walk pieces from v down to the vertex bottom."""
    n = np.zeros(g.nvert, dtype=object)  # exact integers, the counts can be large
    n[bottom] = 1
    lb = int(g.drt[bottom, 0])
    for L in range(lb + 1, g.nlev + 1):
        for v in g.vertices_at_level(L):
            n[v] = sum(n[int(c)] for c in g.down[v] if c >= 0)
    return n


class CommonWalks:
    """The walks that bra and ket share outside a selection of loops, numbered 0..count-1.

    The selection is one complete loop, or one upper and one lower partial loop
    of a crossing loop. Each choice of the common parts embeds the loop in one
    pair of complete walks (bra L, ket R), i.e. one pair of CSFs with the same
    coupling coefficient. The walks are numbered in mixed radix over the
    pieces (the piece nearest the graph head slowest), each piece depth-first
    with the steps d = 0..3 in order.
    """

    def __init__(self, g, loops, mv=None, upper=True, lower=True):
        self.g = g
        self.loops = list(loops)
        segs = []
        for lp in self.loops:
            for sg in common_segments(g, lp, mv, upper, lower):
                if sg not in segs:
                    segs.append(sg)
        self.segments = sorted(segs, key=lambda sg: -int(g.drt[sg[0], 0]))
        self.counts = [walk_counts_to(g, b) for _, b in self.segments]
        self.sizes = [int(n[t]) for (t, _), n in zip(self.segments, self.counts)]
        self.count = 1
        for k in self.sizes:
            self.count *= k

    def arcs(self, k):
        """The arcs (upper vertex, lower vertex, step) of the common walk k (0-based)."""
        if not 0 <= k < self.count:
            raise IndexError(f'common walk {k} out of range [0,{self.count})')
        digits = []
        for size in reversed(self.sizes):
            k, r = divmod(k, size)
            digits.append(r)
        digits.reverse()
        out = []
        for (top, bottom), n, r in zip(self.segments, self.counts, digits):
            v = top
            while v != bottom:
                for d in range(4):
                    c = int(self.g.down[v, d])
                    if c < 0 or n[c] == 0:
                        continue
                    if r < n[c]:
                        out.append((v, c, d))
                        v = c
                        break
                    r -= n[c]
        return out

    def step_vectors(self, k):
        """The complete step vectors (index L-1 for level L) of bra and ket for the common walk k.

        Levels not covered (a half switched off) are -1.
        """
        g = self.g
        bra = np.full(g.nlev, -1, dtype=int)
        ket = np.full(g.nlev, -1, dtype=int)
        for u, _, d in self.arcs(k):
            L = int(g.drt[u, 0])
            bra[L - 1] = ket[L - 1] = d
        for lp in self.loops:
            for n, (db, dk) in enumerate(zip(lp.bra_steps, lp.ket_steps)):
                L = lp.top - n
                bra[L - 1], ket[L - 1] = db, dk
        return bra, ket
