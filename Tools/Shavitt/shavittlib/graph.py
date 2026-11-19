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
"""The split graph (Shavitt graph) of a RASSCF CSF space, rebuilt from the
definition of the CI space in a rasscf.h5 file (or from scratch), and
the decoding of the CSF indices of CI_VECTORS into step vectors and walks.

Conventions (as in src/csf_graph_util/csf_graph.F90):
  * levels 1..nlev; level 0 is the tail; vertex 1 (the head) is at level nlev
  * the step d at level L is the arc from the vertex at level L down to L-1,
    d = 0 (empty), 1 (up, b+1), 2 (down, b-1), 3 (doubly occupied)
  * all indices exposed here are 0-based (vertices, levels-1, CSFs, irreps-1)
    except where noted
"""

from dataclasses import dataclass, field

import numpy as np

STEP_CHAR = '0ud2'


def mul(i, j):
    """Product of two irreps (0-based, D2h and subgroups)."""
    return i ^ j


@dataclass
class Vertex:
    index: int  # 0-based
    level: int
    nel: int
    a: int
    b: int
    c: int
    down: tuple  # 0-based child vertex per step d=0..3, -1: no arc


@dataclass
class ShavittGraph:
    nsym: int
    nactel: int
    ispin: int
    nlev: int
    midlev: int
    mvsta: int  # 0-based first midlevel vertex
    nmidv: int
    lsym: int  # 0-based target irrep of CI_VECTORS
    level_sym: np.ndarray  # irrep (0-based) of each level, index L-1
    level_to_act: np.ndarray  # active orbital (1-based) of each level, index L-1
    drt: np.ndarray  # [nvert,5]: L, N, a, b, c
    down: np.ndarray  # [nvert,4]: 0-based child, -1 none
    now: np.ndarray  # [2,nsym,nmidv]
    iow: np.ndarray  # [2,nsym,nmidv]
    nocsf: np.ndarray  # [nsym(up),nmidv,nsym(tot)]
    iocsf: np.ndarray  # [nsym(up),nmidv,nsym(tot)]
    ncsf: np.ndarray  # [nsym]
    icase: np.ndarray
    npack: int
    nipwlk: int
    nras: np.ndarray = field(default=None)
    nrasel: np.ndarray = field(default=None)
    nrasmx: np.ndarray = field(default=None)
    ci_type: str = ''
    irrep_labels: list = field(default_factory=list)
    ci_vectors: np.ndarray = field(default=None)  # [nroots,nconf]
    root_energies: np.ndarray = field(default=None)
    warnings: list = field(default_factory=list)  # from reading a file

    # ------------------------------------------------------------------
    @classmethod
    def build(cls, nsym, nactel, ispin, lsym, nras, minel=None, maxel=None, level_sym=None, level_to_act=None,
              midlev=None, ci_type=None, irrep_labels=None):
        """The graph of a CSF space from its definition, as SGUGA builds it (see build.py).

        nras[partition][irrep]: orbitals per partition and irrep; minel/maxel:
        cumulative electron bounds of the partitions 1..k (default: CAS);
        lsym: 0-based target irrep. level_sym/level_to_act default to RASSCF's
        numbering (levels partition by partition, irrep by irrep).
        """
        from .build import build_tables, levels_from_partitions
        nras = np.atleast_2d(np.asarray(nras, dtype=int))
        if minel is None:
            if nras.shape[0] != 1:
                raise ValueError('electron bounds needed for more than one partition')
            minel, maxel = [nactel], [-1]
        if level_sym is None or level_to_act is None:
            ls, la = levels_from_partitions(nras)
            level_sym = ls if level_sym is None else level_sym
            level_to_act = la if level_to_act is None else level_to_act
        t = build_tables(nactel, ispin, level_sym, nsym, nras, minel, maxel, midlev)
        self = cls(nsym=nsym, nactel=nactel, ispin=ispin, lsym=lsym, level_sym=np.asarray(level_sym, dtype=int),
                   level_to_act=np.asarray(level_to_act, dtype=int), **t)
        if ci_type is None:
            ci_type = ('CAS' if nras.shape[0] == 1 else
                       'GAS' if maxel is not None and any(int(m) >= 0 for m in maxel) else 'RAS')
        self.ci_type = ci_type
        self.irrep_labels = list(irrep_labels or [])
        return self

    @classmethod
    def from_definition(cls, electrons, spin=1, irrep=1, cas=None, ras1=None, ras2=None, ras3=None, holes=0,
                        particles=0, gas=None, irrep_labels=None):
        """The graph of a CI space given as RASSCF's input defines it (no file).

        electrons: active electrons; spin: 2S+1; irrep: 1-based target irrep.
        One of: cas (orbitals per irrep); ras1/ras2/ras3 (orbitals per irrep)
        with holes (in RAS1) and particles (in RAS3); gas, a list of
        (orbitals per irrep, min, max) with the bounds on the electrons in the
        spaces 1..k together. Raises ValueError for an inconsistent definition.
        """
        from .build import gas_partitions, ras_partitions
        n = int(electrons)
        if n < 0:
            raise ValueError('the number of electrons must be >= 0')
        ras = [list(x) if x else [] for x in (ras1, ras2, ras3)]
        if sum(bool(x) for x in (cas, any(ras), gas)) != 1:
            raise ValueError('give the orbitals of exactly one kind of space: CAS, RAS1..3 or GAS')
        if gas:
            rows = [list(o) for o, _, _ in gas]
        elif any(ras):
            width = max(len(x) for x in ras)
            rows = [x + [0] * (width - len(x)) for x in ras]
        else:
            rows = [list(cas)]
        lens = {len(r) for r in rows}
        if len(lens) != 1 or lens.pop() not in (1, 2, 4, 8):
            raise ValueError('give the orbitals of 1, 2, 4 or 8 irreps, the same number for every space')
        if any(int(x) < 0 for r in rows for x in r):
            raise ValueError('the numbers of orbitals must be >= 0')
        nsym = len(rows[0])
        if not 1 <= irrep <= nsym:
            raise ValueError(f'the target irrep must be 1..{nsym}')
        if gas:
            if int(gas[-1][1]) != n or int(gas[-1][2]) != n:
                raise ValueError(f'GAS: the last space must hold all {n} electrons together (min = max = {n})')
            nras, minel, maxel = gas_partitions(n, rows, [int(g[1]) for g in gas], [int(g[2]) for g in gas])
        elif any(ras):
            if not 0 <= holes <= 2 * sum(rows[0]) or not 0 <= particles <= 2 * sum(rows[2]):
                raise ValueError('RAS: at most 2 holes per RAS1 orbital and 2 particles per RAS3 orbital')
            nras, minel, maxel = ras_partitions(n, rows[0], rows[1], rows[2], int(holes), int(particles))
        else:
            nras, minel, maxel = np.atleast_2d(rows), None, None
        if labels := list(irrep_labels or []):
            if len(labels) != nsym:
                raise ValueError(f'give {nsym} irrep labels')
        return cls.build(nsym, n, int(spin), irrep - 1, nras, minel, maxel, irrep_labels=labels,
                         ci_type='GAS' if gas else None)

    @classmethod
    def from_h5(cls, filename):
        """The graph of the CI space of a rasscf.h5 file, rebuilt from its definition in the root attributes.

        Used: NSYM, NACTEL, SPINMULT, LSYM, NBAS, MO_TYPEINDICES (active
        orbitals), NHOLE1 and NELEC3 (RAS) or GAS_ORBITALS, GAS_MIN_ELECTRONS
        and GAS_MAX_ELECTRONS (GAS), L2ACT (level order) and MIDLEVEL; checked
        against NCONF. Files without MIDLEVEL (or, for GAS, without the
        restrictions) are read with warnings (self.warnings); a graph that
        does not give the CSFs of the file raises ValueError.
        """
        import os

        import h5py
        from .build import gas_partitions, levels_from_partitions, ras_partitions

        if not os.path.isfile(filename):
            raise OSError(f'{filename}: no such file')

        def text(x):
            return x.decode().strip() if isinstance(x, bytes) else str(x).strip()

        warnings = []
        with h5py.File(filename, 'r') as f:
            at = f.attrs
            missing = [k for k in ('NSYM', 'NACTEL', 'SPINMULT', 'LSYM', 'NBAS', 'NCONF') if k not in at]
            if missing or 'MO_TYPEINDICES' not in f:
                raise ValueError(f'{filename}: not a RASSCF wave function file (no {", ".join(missing or ["MO_TYPEINDICES"])})')
            nsym, nactel, ispin = int(at['NSYM']), int(at['NACTEL']), int(at['SPINMULT'])
            lsym = int(at['LSYM']) - 1
            ct = text(at.get('CI_TYPE', b''))
            labels = [text(x) for x in at.get('IRREP_LABELS', [])]
            # active orbitals of each irrep by type (RAS1/2/3; GAS: all 2), MO order irrep by irrep
            types = [text(x).upper() for x in np.asarray(f['MO_TYPEINDICES'][()]).ravel()]
            nbas = np.atleast_1d(np.asarray(at['NBAS'])).astype(int)
            nrs = np.zeros((3, nsym), dtype=int)
            act_sym = []
            i0 = 0
            for s in range(nsym):
                for t in types[i0:i0 + nbas[s]]:
                    if t in ('1', '2', '3'):
                        nrs[int(t) - 1, s] += 1
                        act_sym.append(s)
                i0 += nbas[s]
            act_sym = np.array(act_sym, dtype=int)
            nlev = len(act_sym)
            if nlev == 0:
                raise ValueError(f'{filename}: no active orbitals, no graph')

            if 'GAS_ORBITALS' in f:
                ngssh = np.atleast_2d(np.asarray(f['GAS_ORBITALS'][()])).astype(int)  # [space][irrep]
                nras, minel, maxel = gas_partitions(nactel, ngssh, np.atleast_1d(f['GAS_MIN_ELECTRONS'][()]),
                                                    np.atleast_1d(f['GAS_MAX_ELECTRONS'][()]))
                ci_type = 'GAS'
            elif ct == 'GAS':
                warnings.append('The file has no GAS restrictions (GAS_ORBITALS, GAS_MIN/MAX_ELECTRONS; written by an '
                                'older RASSCF): the graph is built without them, i.e. as a CAS.')
                nras, minel, maxel = nrs.sum(axis=0)[None, :], [nactel], None
                ci_type = 'GAS'
            else:
                nras, minel, maxel = ras_partitions(nactel, nrs[0], nrs[1], nrs[2], int(at.get('NHOLE1', 0)),
                                                    int(at.get('NELEC3', 0)))
                # CI_TYPE says RAS for any CAS with symmetry in older files: use the partitions
                ci_type = 'CAS' if len(nras) == 1 else 'RAS'
            if int(nras.sum()) != nlev:
                raise ValueError(f'{filename}: {int(nras.sum())} orbitals in the CI spaces, but {nlev} active orbitals')

            # the level order: RASSCF numbers the levels partition by partition, irrep by irrep
            level_sym, level_to_act = levels_from_partitions(nras)
            if 'L2ACT' in at:
                l2act = np.asarray(at['L2ACT'][:nlev]).astype(int)
                if sorted(l2act) != list(range(1, nlev + 1)):
                    raise ValueError(f'{filename}: L2ACT is not a permutation of the {nlev} active orbitals')
                if not np.array_equal(l2act, level_to_act):
                    if len(nras) > 1:
                        raise ValueError(f'{filename}: the level order of the file (L2ACT) is not that of the '
                                         'partitions: unknown CSF order')
                    # one partition: any level order is a valid graph
                    level_to_act, level_sym = l2act, act_sym[l2act - 1]
            else:
                warnings.append('The file has no level order (L2ACT): RASSCF\'s default order is assumed.')

            midlev = int(at['MIDLEVEL']) if 'MIDLEVEL' in at else None
            if midlev is None:
                warnings.append('The file has no MIDLEVEL (written by an older RASSCF): it is computed as SGUGA '
                                '(MKMID) does now.')
            self = cls.build(nsym, nactel, ispin, lsym, nras, minel, maxel, level_sym=level_sym,
                             level_to_act=level_to_act, midlev=midlev, ci_type=ci_type, irrep_labels=labels)

            nconf = int(at['NCONF'])
            if 'CI_VECTORS' in f:
                self.ci_vectors = np.atleast_2d(np.asarray(f['CI_VECTORS'][()]))
                nconf_ci = self.ci_vectors.shape[1]
                if nconf_ci != nconf:
                    raise ValueError(f'{filename}: NCONF = {nconf}, but CI_VECTORS has {nconf_ci} coefficients per root')
            if int(self.ncsf[lsym]) != nconf:
                hint = (' (the GAS restrictions are missing from the file)' if ci_type == 'GAS' and 'GAS_ORBITALS' not in f
                        else '')
                raise ValueError(f'{filename}: the rebuilt graph has {int(self.ncsf[lsym])} CSFs of irrep '
                                 f'{self.irrep_name(lsym)}, the file {nconf}{hint}: the CSF order is unknown')
            if 'ROOT_ENERGIES' in f:
                self.root_energies = np.asarray(f['ROOT_ENERGIES'][()])
        self.warnings = warnings
        return self

    # ------------------------------------------------------------------
    @property
    def nvert(self):
        return self.drt.shape[0]

    def vertex(self, iv):
        L, N, a, b, c = (int(x) for x in self.drt[iv])
        return Vertex(iv, L, N, a, b, c, tuple(int(x) for x in self.down[iv]))

    def vertices(self):
        return [self.vertex(iv) for iv in range(self.nvert)]

    def arcs(self):
        """All arcs as (upper vertex, lower vertex, step), 0-based."""
        return [(iv, int(self.down[iv, d]), d) for iv in range(self.nvert) for d in range(4) if self.down[iv, d] >= 0]

    def vertices_at_level(self, L):
        return [iv for iv in range(self.nvert) if self.drt[iv, 0] == L]

    # ------------------------------------------------------------------
    def _unpack(self, pos, levels):
        """Unpack the half walk stored at ICASE[pos:pos+nipwlk] over the given levels (ascending)."""
        steps = []
        word = pos
        ic = int(self.icase[word])
        n = 0
        for _ in levels:
            if n == self.npack:
                n = 0
                word += 1
                ic = int(self.icase[word])
            steps.append(ic & 3)
            ic >>= 2
            n += 1
        return steps

    def step_vector(self, icsf, lsym=None):
        """Step vector (index L-1 for level L) of the 0-based CSF icsf of irrep lsym (0-based)."""
        if lsym is None:
            lsym = self.lsym
        if not 0 <= icsf < self.ncsf[lsym]:
            raise IndexError(f'CSF {icsf} out of range [0,{self.ncsf[lsym]})')
        for mv in range(self.nmidv):
            for isyup in range(self.nsym):
                off = self.iocsf[isyup, mv, lsym]
                n = self.nocsf[isyup, mv, lsym]
                if off <= icsf < off + n:
                    isydwn = mul(lsym, isyup)
                    nup = self.now[0, isyup, mv]
                    loc = icsf - off
                    iup, idwn = loc % nup, loc // nup
                    lower = self._unpack(self.iow[1, isydwn, mv] + idwn * self.nipwlk, range(1, self.midlev + 1))
                    upper = self._unpack(self.iow[0, isyup, mv] + iup * self.nipwlk, range(self.midlev + 1, self.nlev + 1))
                    return np.array(lower + upper, dtype=int)
        raise RuntimeError(f'CSF {icsf} not found in any block')

    def walk(self, steps):
        """Vertices (0-based) visited by a step vector, from the head (level nlev) to the tail (level 0)."""
        path = [0]
        v = 0
        for L in range(self.nlev, 0, -1):
            v = int(self.down[v, steps[L - 1]])
            if v < 0:
                raise ValueError(f'step vector not a walk of the graph (level {L})')
            path.append(v)
        return path

    def csf_string(self, steps, by_active=False):
        """The step vector as 0/u/d/2, in level order or in active-orbital order."""
        s = [STEP_CHAR[d] for d in steps]
        if by_active:
            out = [''] * self.nlev
            for L in range(self.nlev):
                out[self.level_to_act[L] - 1] = s[L]
            s = out
        return ''.join(s)

    def leading_csfs(self, root=0, threshold=0.05):
        """(index, coefficient, step vector) of the CSFs of a root with |c| >= threshold, largest first."""
        if self.ci_vectors is None:
            raise ValueError('no CI_VECTORS in the file')
        c = self.ci_vectors[root]
        idx = np.nonzero(np.abs(c) >= threshold)[0]
        idx = idx[np.argsort(-np.abs(c[idx]))]
        return [(int(i), float(c[i]), self.step_vector(int(i))) for i in idx]

    # ------------------------------------------------------------------
    @staticmethod
    def parse_steps(text):
        """A step vector from a string of 0/u/d/2 (or 0/1/2/3 as steps with a leading 'd:'), blanks ignored."""
        t = ''.join(text.split())
        if t.lower().startswith('d:'):
            return np.array([int(ch) for ch in t[2:]], dtype=int)
        try:
            return np.array([STEP_CHAR.index(ch) for ch in t.lower()], dtype=int)
        except ValueError:
            raise ValueError(f"'{text}': use the characters 0, u, d, 2") from None

    def walk_irrep(self, steps, levels=None):
        """Irrep (0-based) of a (partial) step vector over the given levels (1-based, default all)."""
        if levels is None:
            levels = range(1, self.nlev + 1)
        s = 0
        for L in levels:
            if steps[L - 1] in (1, 2):
                s = mul(s, int(self.level_sym[L - 1]))
        return s

    def csf_index(self, steps):
        """0-based CSF index of a complete walk in CI_VECTORS, or None if it is not of the target irrep."""
        steps = np.asarray(steps, dtype=int)
        if len(steps) != self.nlev:
            raise ValueError(f'a step vector has {self.nlev} steps')
        path = self.walk(steps)
        if self.walk_irrep(steps) != self.lsym:
            return None
        mv = path[self.nlev - self.midlev] - self.mvsta
        isyup = self.walk_irrep(steps, range(self.midlev + 1, self.nlev + 1))
        isydwn = mul(self.lsym, isyup)
        nup = self.now[0, isyup, mv]
        ndwn = self.now[1, isydwn, mv]
        lo = list(steps[:self.midlev])
        up = list(steps[self.midlev:])
        idwn = next(i for i in range(ndwn)
                    if self._unpack(self.iow[1, isydwn, mv] + i * self.nipwlk, range(1, self.midlev + 1)) == lo)
        iup = next(i for i in range(nup)
                   if self._unpack(self.iow[0, isyup, mv] + i * self.nipwlk, range(self.midlev + 1, self.nlev + 1)) == up)
        return int(self.iocsf[isyup, mv, self.lsym] + idwn * nup + iup)

    def midvertices(self):
        """0-based indices of the vertices at the midlevel."""
        return self.vertices_at_level(self.midlev)

    def path_counts(self):
        """(number of walks head->v, number of walks v->tail) for all vertices, all irreps together."""
        order = sorted(range(self.nvert), key=lambda v: -self.drt[v, 0])
        up = np.zeros(self.nvert, dtype=np.int64)
        up[0] = 1
        for v in order:
            for d in range(4):
                c = self.down[v, d]
                if c >= 0:
                    up[c] += up[v]
        dn = np.zeros(self.nvert, dtype=np.int64)
        for v in reversed(order):
            if self.drt[v, 0] == 0:
                dn[v] = 1
            for d in range(4):
                c = self.down[v, d]
                if c >= 0:
                    dn[v] += dn[c]
        return up, dn

    def arcs_through(self, v):
        """The arcs (upper, lower, step) on the walks that pass through vertex v."""
        below = {v}
        for L in range(self.drt[v, 0], 0, -1):
            for u in [u for u in below if self.drt[u, 0] == L]:
                below.update(int(c) for c in self.down[u] if c >= 0)
        above = {v}
        for L in range(self.drt[v, 0] + 1, self.nlev + 1):
            for u in self.vertices_at_level(L):
                if any(c in above for c in self.down[u] if c >= 0):
                    above.add(u)
        return [(u, c, d) for u, c, d in self.arcs() if (u in below and c in below) or (u in above and c in above)]

    # ------------------------------------------------------------------
    def count_walks(self):
        """Number of complete walks by total irrep, enumerated directly from the graph (cross-check of NCSF)."""
        # weights[v][s]: number of walks from v down to the tail with accumulated irrep s
        w = np.zeros((self.nvert, self.nsym), dtype=np.int64)
        order = sorted(range(self.nvert), key=lambda v: self.drt[v, 0])
        for v in order:
            if self.drt[v, 0] == 0:
                w[v, 0] = 1
                continue
            L = self.drt[v, 0]
            for d in range(4):
                c = self.down[v, d]
                if c < 0:
                    continue
                s = self.level_sym[L - 1] if d in (1, 2) else 0
                for t in range(self.nsym):
                    w[v, mul(s, t)] += w[c, t]
        return w[0]

    def irrep_name(self, isym):
        """Label of the 0-based irrep isym."""
        if isym < len(self.irrep_labels) and self.irrep_labels[isym]:
            return self.irrep_labels[isym]
        return str(isym + 1)

    def csf_counts(self, sep=', '):
        """The number of CSFs of each irrep, e.g. 'a1 139 (wave function), b1 120, ...'."""
        return sep.join(f'{self.irrep_name(i)} {int(n)}' + (' (wave function)' if i == self.lsym else '')
                        for i, n in enumerate(self.ncsf))

    def summary(self):
        lines = [
            f'CI type {self.ci_type}, {self.nactel} electrons in {self.nlev} levels, 2S+1 = {self.ispin}, '
            f'{self.nsym} irreps, target irrep {self.irrep_name(self.lsym)}',
            f'{self.nvert} vertices, {len(self.arcs())} arcs, midlevel {self.midlev} ({self.nmidv} midvertices)',
            'CSFs per irrep: ' + self.csf_counts(),
            'level  irrep  act   vertices',
        ]
        for L in range(self.nlev, 0, -1):
            lines.append(f'{L:5d}  {self.irrep_name(self.level_sym[L - 1]):>5s}  {self.level_to_act[L - 1]:3d}   '
                         f'{len(self.vertices_at_level(L))}')
        return '\n'.join(lines)
