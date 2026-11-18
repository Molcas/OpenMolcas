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
definition in the SHAVITT group of a rasscf.h5 file (or from scratch), and
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
    def from_h5(cls, filename):
        import h5py
        with h5py.File(filename, 'r') as f:
            if 'SHAVITT' not in f:
                raise KeyError(f'{filename}: no SHAVITT group; rerun RASSCF with a version that writes the split graph')
            g = f['SHAVITT']
            at = g.attrs

            # HDF5 stores the Fortran arrays with the dimensions reversed
            def F(name):
                return np.asarray(g[name][()]).T

            ct = f.attrs.get('CI_TYPE', b'')
            ct = ct.decode() if isinstance(ct, bytes) else str(ct)
            labels = [x.decode().strip() if isinstance(x, bytes) else str(x).strip()
                      for x in f.attrs.get('IRREP_LABELS', [])]
            nsym = int(f.attrs['NSYM'])
            lsym = int(f.attrs['LSYM']) - 1
            if 'DRT' not in g:
                # the definition of the graph: rebuild it as SGUGA does, check it against RASSCF's counts
                nras = F('PARTITION_ORBITALS').astype(int).T  # [partition][irrep]
                level_sym = np.atleast_1d(F('LEVEL_IRREPS')).astype(int) - 1
                nlev = len(level_sym)
                self = cls.build(nsym, int(f.attrs['NACTEL']), int(f.attrs['SPINMULT']), lsym, nras,
                                 np.atleast_1d(F('PARTITION_MIN_ELECTRONS')).astype(int),
                                 np.atleast_1d(F('PARTITION_MAX_ELECTRONS')).astype(int), level_sym=level_sym,
                                 level_to_act=np.asarray(f.attrs['L2ACT'][:nlev]).astype(int),
                                 midlev=int(at['MIDLEVEL']), irrep_labels=labels)
                ncsf = np.atleast_1d(F('CSFS_PER_IRREP')).astype(int)
                if self.nvert != int(at['NUM_VERTICES']) or not np.array_equal(self.ncsf, ncsf):
                    raise ValueError(f'{filename}: the rebuilt graph ({self.nvert} vertices, '
                                     f'CSFs {self.ncsf.tolist()}) '
                                     f'differs from RASSCF\'s ({int(at["NUM_VERTICES"])} vertices, CSFs '
                                     f'{ncsf.tolist()}): the CSF order of SGUGA may have changed')
            else:
                # files written by Shavitt 0.9 (before 2026-11-18): the tables as stored by RASSCF
                down = F('DOWN').astype(int) - 1
                self = cls(
                    nsym=nsym, nactel=int(at['NACTEL']), ispin=int(at['ISPIN']), nlev=int(at['NLEV']),
                    midlev=int(at['MIDLEV']), mvsta=int(at['MVSTA']) - 1, nmidv=int(at['NMIDV']), lsym=lsym,
                    level_sym=F('LEVEL_SYM').astype(int) - 1, level_to_act=F('LEVEL_TO_ACT').astype(int),
                    drt=F('DRT').astype(int), down=down,
                    now=F('NOW').astype(int), iow=F('IOW').astype(int),
                    nocsf=F('NOCSF').astype(int), iocsf=F('IOCSF').astype(int), ncsf=F('NCSF').astype(int),
                    icase=F('ICASE').astype(np.int64), npack=int(at['NPACK']), nipwlk=int(at['NIPWLK']),
                    nras=F('NRAS').astype(int), nrasel=F('NRASEL').astype(int), nrasmx=F('NRASMX').astype(int),
                )
                # CI_TYPE says RAS for any CAS with symmetry in files before that fix: use the partitions
                self.ci_type = 'CAS' if int(at['NRSPRT']) == 1 else 'RAS'
                if 'L2ACT' in f.attrs:  # level -> active orbital (files written before LEVEL_TO_ACT was IDXSX)
                    self.level_to_act = np.asarray(f.attrs['L2ACT'][:self.nlev]).astype(int)
                self.irrep_labels = labels
            if ct.strip() == 'GAS':
                self.ci_type = 'GAS'
            if 'CI_VECTORS' in f:
                self.ci_vectors = np.atleast_2d(np.asarray(f['CI_VECTORS'][()]))
            if 'ROOT_ENERGIES' in f:
                self.root_energies = np.asarray(f['ROOT_ENERGIES'][()])
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

    def summary(self):
        lines = [
            f'CI type {self.ci_type}, {self.nactel} electrons in {self.nlev} levels, 2S+1 = {self.ispin}, '
            f'{self.nsym} irreps, target irrep {self.irrep_name(self.lsym)}',
            f'{self.nvert} vertices, {len(self.arcs())} arcs, midlevel {self.midlev} ({self.nmidv} midvertices)',
            'CSFs per irrep: ' + ' '.join(str(int(n)) for n in self.ncsf),
            'level  irrep  act   vertices',
        ]
        for L in range(self.nlev, 0, -1):
            lines.append(f'{L:5d}  {self.irrep_name(self.level_sym[L - 1]):>5s}  {self.level_to_act[L - 1]:3d}   '
                         f'{len(self.vertices_at_level(L))}')
        return '\n'.join(lines)
