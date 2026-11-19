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
"""The segments of the singlet generators at a vertex, as in SGUGA (MKSEG).

A segment spans one level: a top pair of vertices (left/bra, right/ket) at
level L, the bra and ket steps d at that level, and the bottom pair at L-1.
SGUGA tabulates them by the top-left vertex (ISGM(IVLT,ISGT), VSGM): the
right top vertex is the same vertex (top class 0: above a loop, 3: below) or
its partner of the same level (class 1: a_R = a_L, b_R = b_L - 1; class 2:
a_R = a_L - 1, b_R = b_L + 1). The value is a formula (ISVC) in b of the
top-left vertex (ASTA segment values).

Catalog (SGUGA numbering 1..34):
  1- 4  walk segments above a loop       23-26  walk segments below a loop
  5- 8  head segments (level j of E_ij)  19-22  tail segments (level i)
  9-13  intermediate, db = -1            14-18  intermediate, db = +1
 27-30  weight segments of E_ii above    31-34  weight segments of E_ii below
"""

from dataclasses import dataclass

from .loops import pair_class, segment_value

# segments 1..34 of sguga.F90 (index 0..33 here)
ITVPT = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 3, 3, 3, 3]
IBVPT = [0, 0, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 2, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 3, 3, 3, 3]
IC1 = [0, 1, 2, 3, 0, 2, 0, 1, 0, 1, 1, 2, 3, 0, 1, 2, 2, 3, 1, 3, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
IC2 = [0, 1, 2, 3, 1, 3, 2, 3, 0, 1, 2, 2, 3, 0, 1, 1, 2, 3, 0, 2, 0, 1, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
ISVC = [1, 1, 1, 1, 1, 7, 8, 4, 1, 2, 9, 10, 2, 1, 2, 11, 12, 2, 1, 5, 6, 3, 1, 1, 1, 1, 0, 1, 1, 13, 0, 1, 1, 13]
NSEG = len(ISVC)

FORMULA = {0: '0', 1: '1', 2: '-1', 3: '(-1)^b', 4: '(-1)^(b+1)', 5: 'sqrt(b/(b+1))',
           6: '(-1)^b sqrt((b+1)/(b+2))', 7: 'sqrt((b+2)/(b+1))', 8: '(-1)^b sqrt((b+2)/(b+1))',
           9: '(-1)^b / sqrt(b(b+1))', 10: 'sqrt(b(b+2))/(b+1)', 11: '(-1)^b / sqrt((b+1)(b+2))',
           12: 'sqrt((b+1)(b+3))/(b+2)', 13: '2'}

CLASS = {0: 'same vertex (above the loop)', 1: 'partner, db = -1 (a_R = a_L, b_R = b_L - 1)',
         2: 'partner, db = +1 (a_R = a_L - 1, b_R = b_L + 1)', 3: 'same vertex (below the loop)'}

# category of each segment: (key, label)
CATEGORY = [('walk', 'walk above')] * 4 + [('head', 'head')] * 4 + [('mid', 'mid db=-1')] * 5 + \
           [('mid', 'mid db=+1')] * 5 + [('tail', 'tail')] * 4 + [('walk', 'walk below')] * 4 + \
           [('weight', 'weight above')] * 4 + [('weight', 'weight below')] * 4

DESCRIPTION = {
    'walk above': 'bra and ket on the same walk above the loop head; value 1',
    'head': 'the loop head: level j of E_ij, where bra and ket separate',
    'mid db=-1': 'inside the loop, between head and tail; bra and ket differ by db = -1',
    'mid db=+1': 'inside the loop, between head and tail; bra and ket differ by db = +1',
    'tail': 'the loop tail: level i of E_ij, where bra and ket join again',
    'walk below': 'bra and ket on the same walk below the loop tail; value 1',
    'weight above': 'diagonal operator E_ii (occupation of the level), upper walk',
    'weight below': 'diagonal operator E_ii (occupation of the level), lower walk',
}


@dataclass
class Segment:
    k: int  # 0-based index in the catalog (SGUGA number k+1)
    vlt: int  # top left (bra) vertex
    vrt: int  # top right (ket) vertex
    vlb: int  # bottom left
    vrb: int  # bottom right
    value: float
    b: int  # b of the top-left vertex, the argument of the formula
    role: str  # 'L': the vertex asked for is the top-left vertex, 'R': the top-right one

    @property
    def number(self):
        return self.k + 1

    @property
    def category(self):
        return CATEGORY[self.k][1]

    @property
    def kind(self):
        return CATEGORY[self.k][0]

    @property
    def steps(self):
        return '0ud2'[IC1[self.k]], '0ud2'[IC2[self.k]]

    @property
    def formula(self):
        return FORMULA[ISVC[self.k]]


def partner(g, v, c):
    """The vertex of the same level that forms a pair of class c (1, 2) with v as the left vertex, or None."""
    for u in g.vertices_at_level(int(g.drt[v, 0])):
        if pair_class(g, v, u) == c:
            return u
    return None


def _segment(g, k, vl, vr, role):
    lb, rb = int(g.down[vl, IC1[k]]), int(g.down[vr, IC2[k]])
    if lb < 0 or rb < 0:
        return None
    b = int(g.drt[vl, 3])
    return Segment(k, vl, vr, lb, rb, float(segment_value(ISVC[k], b)), b, role)


def segments_at(g, v, as_left=True, as_right=True, kinds=None):
    """The segments with v as top-left (bra) and/or top-right (ket) vertex, in SGUGA order.

    kinds: the categories to keep ('walk', 'head', 'mid', 'tail', 'weight'), default all.
    """
    out = []
    if int(g.drt[v, 0]) == 0:
        return out  # the tail of the graph has no segments below it
    if as_left:
        for k in range(NSEG):
            t = ITVPT[k]
            vr = v if t in (0, 3) else partner(g, v, t)
            if vr is None:
                continue
            s = _segment(g, k, v, vr, 'L')
            if s is not None:
                out.append(s)
    if as_right:
        # v as the ket vertex of a pair of class 1 or 2 (classes 0 and 3 are the same vertex: listed above)
        for u in g.vertices_at_level(int(g.drt[v, 0])):
            c = pair_class(g, u, v)
            if c not in (1, 2):
                continue
            for k in range(NSEG):
                if ITVPT[k] == c:
                    s = _segment(g, k, u, v, 'R')
                    if s is not None:
                        out.append(s)
    if kinds is not None:
        out = [s for s in out if s.kind in kinds]
    return out
