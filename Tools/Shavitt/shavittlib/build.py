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
"""Construction of the split graph (Shavitt graph) of a CAS/RAS/GAS CSF space
from its definition, as SGUGA does it (src/csf_graph_util/csf_graph.F90,
MKSGUGA), so that the CSFs come out in the order of CI_VECTORS:

  mkDRT0  the full Paldus distinct-row table, built level by level from the
          head (a0, b0, c0), new vertices numbered in the order they are met
          (upper vertex, then step d = 0..3)
  RMVERT  RAS/GAS: remove the vertices that violate the cumulative electron
          bounds at the partition boundaries, and those no longer on a walk
          from the head to the tail; renumber keeping the order
  MKMID   the midlevel: the level with the most balanced number of upper and
          lower walks (first such level from the bottom)
  MKCOT   the upper and lower half walks, depth first with the steps 0..3
          tried in order; counted and numbered by half, irrep and midvertex
  CSFCOUNT  CSF blocks: midvertex outermost, upper-walk irrep inner; within a
          block the upper half walk runs fastest

Levels are numbered partition by partition, within a partition irrep by
irrep (as RASSCF does); active orbitals irrep by irrep, within an irrep
partition by partition.
"""

import numpy as np

# change of (a, b, c) going down one level with the step d (Steps of mkDRT0)
DSTEP = ((0, 0, 1), (0, 1, 0), (1, -1, 1), (1, 0, 0))
NPACK = 31  # steps per ICASE word, as SGUGA with 64-bit integers


def _tdiv(a, b):
    """Integer division truncated towards zero (Fortran)."""
    q = abs(a) // abs(b)
    return q if (a >= 0) == (b > 0) else -q


def levels_from_partitions(nras):
    """Level irreps (0-based) and the level -> active orbital map (1-based) for nras[partition][irrep]."""
    nras = np.asarray(nras, dtype=int)
    npart, nsym = nras.shape
    # active orbitals: irrep by irrep, within an irrep partition by partition
    act = {}
    n = 0
    for s in range(nsym):
        for p in range(npart):
            for k in range(nras[p, s]):
                n += 1
                act[(p, s, k)] = n
    level_sym, level_to_act = [], []
    for p in range(npart):
        for s in range(nsym):
            for k in range(nras[p, s]):
                level_sym.append(s)
                level_to_act.append(act[(p, s, k)])
    return np.array(level_sym, dtype=int), np.array(level_to_act, dtype=int)


def ras_partitions(nactel, nrs1, nrs2, nrs3, nhole1, nelec3):
    """Partitions and cumulative bounds of a RAS (CAS) space, as Setup_RASSCF.

    nrs1..3: orbitals of each irrep; returns nras[partition][irrep], minel,
    maxel (None: no upper bounds). Without holes and particles one partition.
    """
    nrs1, nrs2, nrs3 = (np.asarray(x, dtype=int) for x in (nrs1, nrs2, nrs3))
    if nhole1 + nelec3 == 0:
        return np.atleast_2d(nrs1 + nrs2 + nrs3), [nactel], None
    return np.array([nrs1, nrs2, nrs3]), [2 * int(nrs1.sum()) - nhole1, nactel - nelec3, nactel], None


def gas_partitions(nactel, ngssh, gmin, gmax):
    """Partitions and cumulative bounds of a GAS space, as Space_From_GAS (src/csf_graph_util/ci_space_def.F90).

    ngssh[space][irrep]: orbitals; gmin/gmax: bounds on the electrons in the
    spaces 1..g together. Spaces without orbitals are dropped (their bounds
    restrict the previous partition), the bounds are tightened to the possible
    ones and the last partition holds all electrons.
    """
    ngssh = np.atleast_2d(np.asarray(ngssh, dtype=int))
    nras, mn, mx = [], [], []
    for g, row in enumerate(ngssh):
        if row.sum() == 0:
            if not nras:
                if gmin[g] > 0 or gmax[g] < 0:
                    raise ValueError(f'GAS space {g + 1}: inconsistent electron bounds')
            else:
                mn[-1] = max(mn[-1], int(gmin[g]))
                mx[-1] = min(mx[-1], int(gmax[g]))
            continue
        nras.append(row)
        mn.append(int(gmin[g]))
        mx.append(int(gmax[g]))
    if not nras:
        raise ValueError('GAS: no active orbitals')
    nras = np.array(nras)
    sizes = nras.sum(axis=1)
    for k in range(len(nras)):
        below, above = int(sizes[:k + 1].sum()), int(sizes[k + 1:].sum())
        mn[k] = max(mn[k], 0, nactel - 2 * above)
        mx[k] = min(mx[k], nactel, 2 * below)
    mn[-1] = mx[-1] = nactel
    for k in range(len(nras)):
        if mn[k] > mx[k]:
            raise ValueError(f'GAS partition {k + 1}: inconsistent electron bounds ({mn[k]} > {mx[k]})')
    return nras, mn, mx


def build_tables(nactel, ispin, level_sym, nsym, nras=None, minel=None, maxel=None, midlev=None):
    """The graph tables of SGUGA for nactel electrons, spin multiplicity ispin and the level irreps.

    nras[partition][irrep]: orbitals per partition and irrep (levels numbered
    partition by partition); minel/maxel: lower/upper bound on the electrons in
    the partitions 1..k (cumulative, maxel -1: no bound). None: one partition.
    midlev: force the midlevel (None: MKMID).
    Returns a dict with the 0-based arrays of ShavittGraph.
    """
    level_sym = np.asarray(level_sym, dtype=int)
    nlev = len(level_sym)
    b0 = ispin - 1
    a0 = (nactel - b0) // 2
    c0 = nlev - a0 - b0
    if 2 * a0 + b0 != nactel or min(a0, b0, c0) < 0:
        raise ValueError(f'impossible graph: {nactel} electrons, {nlev} orbitals, 2S+1 = {ispin}')
    if nlev < 1:
        raise ValueError('no active orbitals: no graph')

    # mkDRT0: full Paldus table, vertices in order of creation
    abc = [(a0, b0, c0)]
    index = {abc[0]: 0}
    down = [[-1] * 4]
    vsta, vend = 0, 1
    for _ in range(nlev, 0, -1):
        for vu in range(vsta, vend):
            a, b, c = abc[vu]
            for d, (da, db, dc) in enumerate(DSTEP):
                t = (a - da, b - db, c - dc)
                if min(t) < 0:
                    continue
                if t not in index:
                    index[t] = len(abc)
                    abc.append(t)
                    down.append([-1] * 4)
                down[vu][d] = index[t]
        vsta, vend = vend, len(abc)
    nv0 = len(abc)
    lev = [a + b + c for a, b, c in abc]
    nel = [2 * a + b for a, b, c in abc]

    # RMVERT: electron bounds at the partition boundaries
    if nras is None:
        nras = [np.bincount(level_sym, minlength=nsym)]
        minel, maxel = [nactel], [-1]
    nras = np.asarray(nras, dtype=int)
    lim = [0] * (nlev + 1)
    limmx = [10**9] * (nlev + 1)
    L = 0
    for p in range(nras.shape[0]):
        L += int(nras[p].sum())
        if L > 0:
            lim[L] = int(minel[p])
            if maxel is not None and int(maxel[p]) >= 0:
                limmx[L] = int(maxel[p])
    ver = [True] * nv0
    for iv in range(nv0 - 1):
        if nel[iv] < lim[lev[iv]] or nel[iv] > limmx[lev[iv]]:
            ver[iv] = False
    while True:
        changes = 0
        for iv in range(nv0 - 1):
            if not ver[iv]:
                down[iv] = [-1] * 4
                continue
            for d in range(4):
                c = down[iv][d]
                if c >= 0 and not ver[c]:
                    down[iv][d] = -1
                    changes += 1
        for iv in range(nv0 - 1):
            if ver[iv] and all(c < 0 for c in down[iv]):
                ver[iv] = False
                changes += 1
        conn = [False] * nv0
        conn[0] = ver[0]
        for iv in range(nv0 - 1):
            if ver[iv]:
                for c in down[iv]:
                    if c >= 0 and ver[c]:
                        conn[c] = True
        for iv in range(nv0):
            if ver[iv] and not conn[iv]:
                ver[iv] = False
                changes += 1
        if changes == 0:
            break
    if not ver[0]:
        raise ValueError('too severe restrictions: no CSFs')
    new = {}
    for iv in range(nv0):
        if ver[iv]:
            new[iv] = len(new)
    nvert = len(new)
    drt = np.zeros((nvert, 5), dtype=int)
    dn = np.full((nvert, 4), -1, dtype=int)
    for iv, jv in new.items():
        a, b, c = abc[iv]
        drt[jv] = (lev[iv], nel[iv], a, b, c)
        for d in range(4):
            c = down[iv][d]
            if c >= 0:
                dn[jv, d] = new[c]

    # MKDAW/MKRAW: walks below (to the tail) and above (from the head) each vertex
    daw = np.zeros(nvert, dtype=object)
    daw[nvert - 1] = 1
    for iv in range(nvert - 2, -1, -1):
        daw[iv] = sum(daw[c] for c in dn[iv] if c >= 0)
    raw = np.zeros(nvert, dtype=object)
    raw[0] = 1
    for iv in range(nvert):
        for c in dn[iv]:
            if c >= 0:
                raw[c] += raw[iv]

    # MKMID
    if midlev is None:
        midlev = min(nlev, 1)
        minw = 1000000
        for il in range(1, nlev):
            nw = sum(abs(raw[iv] - daw[iv]) for iv in range(nvert) if drt[iv, 0] == il)
            if nw < minw:
                midlev, minw = il, nw
    mids = [iv for iv in range(nvert) if drt[iv, 0] == midlev]
    mvsta, nmidv = mids[0], len(mids)

    # MKCOT: half walks depth first, steps 0..3
    nipwlk = max(1 + _tdiv(midlev - 1, NPACK), 1 + _tdiv(nlev - midlev - 1, NPACK))
    walks = {}  # (half, irrep, mv) -> list of step lists (levels ascending)

    def dfs(v, L, bottom, sym, steps, half):
        if L == bottom:
            mv = (v if half == 1 else start) - mvsta
            walks.setdefault((half, sym, mv), []).append(list(reversed(steps)))
            return
        for d in range(4):
            c = dn[v, d]
            if c < 0:
                continue
            s = mul(sym, level_sym[L - 1]) if d in (1, 2) else sym
            dfs(int(c), L - 1, bottom, s, steps + [d], half)

    start = 0
    dfs(0, nlev, midlev, 0, [], 1)
    for start in mids:
        dfs(start, midlev, 0, 0, [], 2)
    now = np.zeros((2, nsym, nmidv), dtype=int)
    for (h, s, mv), ws in walks.items():
        now[h - 1, s, mv] = len(ws)
    iow = np.zeros((2, nsym, nmidv), dtype=int)
    nw = 0
    for h in range(2):
        for mv in range(nmidv):
            for s in range(nsym):
                iow[h, s, mv] = nw * nipwlk
                nw += now[h, s, mv]
    icase = np.zeros(nw * nipwlk, dtype=np.int64)
    for (h, s, mv), ws in walks.items():
        for k, steps in enumerate(ws):
            pos = iow[h - 1, s, mv] + k * nipwlk
            for w in range(0, len(steps), NPACK):
                ic = 0
                for d in reversed(steps[w:w + NPACK]):
                    ic = 4 * ic + d
                icase[pos + w // NPACK] = ic

    # CSFCOUNT
    nocsf = np.zeros((nsym, nmidv, nsym), dtype=int)
    iocsf = np.zeros((nsym, nmidv, nsym), dtype=int)
    ncsf = np.zeros(nsym, dtype=int)
    for st in range(nsym):
        for mv in range(nmidv):
            for su in range(nsym):
                n = now[0, su, mv] * now[1, mul(st, su), mv]
                nocsf[su, mv, st] = n
                iocsf[su, mv, st] = ncsf[st]
                ncsf[st] += n

    return dict(nlev=nlev, midlev=midlev, mvsta=mvsta, nmidv=nmidv, drt=drt, down=dn, now=now, iow=iow,
                nocsf=nocsf, iocsf=iocsf, ncsf=ncsf, icase=icase, npack=NPACK, nipwlk=nipwlk,
                nras=nras.T.copy(), nrasel=np.asarray(minel, dtype=int),
                nrasmx=np.asarray(maxel if maxel is not None else [-1] * nras.shape[0], dtype=int))


def mul(i, j):
    """Product of two irreps (0-based, D2h and subgroups)."""
    return i ^ j
