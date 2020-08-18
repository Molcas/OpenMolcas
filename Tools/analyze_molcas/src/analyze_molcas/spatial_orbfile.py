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
# Copyright (C) 2020, Oskar Weser                                      *
#***********************************************************************

from typing import Sequence, TextIO
from os import PathLike
from copy import deepcopy

import numpy as np
from numpy import argsort, isclose, reshape, argwhere
from scipy.linalg import block_diag

from analyze_molcas._base_orbfile import (
    _Orbitals, FileFormat, _reindex,
    _read_spatial_orbs, _read_header,
    _write_section_1D, _write_section_2D)

from analyze_molcas.spin_orbfile import SpinOrbs

from analyze_molcas._linalg_helper import get_kabsch_transformation


class SpatialOrbs(_Orbitals):
    """Represents a RasOrb file with spatial orbitals

The operators +, *, ==, and != are overloaded.
Similar to lists `+` concatenates two SpatialOrbs and
`*` concatenates a SpatialOrb with itself n-times.


Attributes
----------
All attributes are lists, whose length is the number of
irreps.  The i-th element is the i-th irrep.

orbs : list of integer
    The number of orbitals per irrep.
    `RasOrb.orb[i]` is the number of orbitals in the i-th irrep.
    `sum(RasOrb.orb)` is the total number of orbitals.
coeff : list of float 2D arrays
    The MO-coefficients per irrep.
    `RasOrb.coeff[i][mu, l]` is the linear coefficient of the
    mu-th AO in the l-th MO in the i-th irrep.
occ : list of float
    The occupation numbers per irrep.
    `RasOrb.occ[i][j]` is the occupation number of the j-th MO
    in the i-th irrep.
energy : list of float
    The orbital energies per irrep.
    `RasOrb.energy[i][j]` is the energy of the j-th MO
    in the i-th irrep.
idx : list of string 1D arrays
    The RAS index per irrep.
    `RasOrb.energy[i][j]` is the RAS index of the j-th MO
    in the i-th irrep.

Examples
--------
Read a RasOrb file:

>>> ras_orb = SpatialOrbs.read_orbfile('./quint.RasOrb')

Write a RasOrb file:

>>> ras_orb.write_orbfile('./transformed.RasOrb')

Put all orbitals in active space:

>>> for index in ras_orb.idx:
>>>     index[:] = '2'
    """

    @classmethod
    def read_orbfile(cls, path: PathLike):
        with open(path, 'r') as f:
            return cls.read_orbfile_stream(f)

    @classmethod
    def read_orbfile_stream(cls, stream: TextIO):
        version, orbs, spin_orbs = _read_header(stream)

        if spin_orbs:
            raise FileFormat('File represents spin orbitals')

        coeff, occ, energy, idx = _read_spatial_orbs(stream, orbs, version)
        return cls(orbs, coeff, occ, energy, idx)

    def reindex(self,
                new_idx: Sequence[Sequence[int]], inplace: bool=False):
        if inplace:
            self.coeff = [
                coeff[:, idx] for idx, coeff in zip(new_idx, self.coeff)]
            self.energy = _reindex(self.energy, new_idx)
            self.idx = _reindex(self.idx, new_idx)
            self.occ = _reindex(self.occ, new_idx)
        else:
            new = self.copy()
            new.reindex(new_idx, inplace=True)
            return new

    def round_occ(self, inplace: bool=False):
        new_occ = [occ.round() for occ in self.occ]

        if not isclose(sum(occ.sum() for occ in self.occ),
                       sum(occ.sum() for occ in new_occ)):
            raise ValueError('Rounded occupation numbers yield different '
                             'number of electrons.')
        if inplace:
            self.occ = new_occ
        else:
            new = self.copy()
            new.round_occ(inplace=True)
            return new

    def to_SpinOrbs(self) -> SpinOrbs:
        spat = self.reindex([argsort(-v, kind='stable') for v in self.occ])

        new_occ = {'a': [occ / 2.0 for occ in spat.occ],
                   'b': [occ / 2.0 for occ in spat.occ]}

        return SpinOrbs(
            orbs=spat.orbs.copy(), coeff=self._spin_copy(spat.coeff),
            occ=new_occ, energy=self._spin_copy(spat.energy),
            idx=spat.idx.copy())

    @staticmethod
    def _spin_copy(v):
        return {'a': deepcopy(v), 'b': deepcopy(v)}

    def _write_header(self):
        return ['#INPORB 2.2',
                '#INFO',
                '* PyOrb written',
                f"{0:8d}{len(self.orbs):8d}{0:8d}",
                f"{''.join(f'{x:8d}' for x in self.orbs)}",
                f"{''.join(f'{x:8d}' for x in self.orbs)}",
                '* SpatialOrbs written from analyze_molcas']

    def _write_MO_coeff(self):
        def get_paragraph(irrep, orb):
            return f'* ORBITAL{irrep + 1:5d}{orb + 1:5d}'

        yield from _write_section_2D(
            title='#ORB', paragraph=get_paragraph,
            orbs=self.orbs, values=self.coeff)

    def _write_MO_occ(self):
        yield from _write_section_1D(
            title='#OCC', subtitle='* OCCUPATION NUMBERS',
            orbs=self.orbs, values=self.occ)

    def _write_MO_human_occ(self):
        yield from _write_section_1D(
            title='#OCHR', subtitle='* OCCUPATION NUMBERS (HUMAN-READABLE)',
            orbs=self.orbs, values=self.occ, cols=10, fmt='8.4f')

    def _write_MO_E(self):
        yield from _write_section_1D(
            title='#ONE', subtitle='* ONE ELECTRON ENERGIES',
            orbs=self.orbs, values=self.energy, cols=10, fmt='12.4E')

    def canonicalize(self):
        """Sort into canonical MOLCAS order."""
        def to_key(v):
            ORB_INDEX = ['f', 'i', '1', '2', '3', 's', 'd']
            D = {index: i for i, index in enumerate(ORB_INDEX)}
            return np.array([D[i] for i in v], dtype='i8')

        idx = [argsort(to_key(self.idx[0]), kind='stable')]
        new = self.reindex(idx)
        return new

    def get_idx_active_space(self):
        return [reshape(argwhere(condition), sum(condition))
                for condition in (index == '2' for index in self.idx)]

    def __add__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError("Can only add {}".format(self.__class__))
        if len(self.orbs) != len(other.orbs):
            raise ValueError("Symmetry group does not match.")
        return self.__class__(
            orbs=[n_s + n_o for n_s, n_o in zip(self.orbs, other.orbs)],
            coeff=[block_diag(Q, P) for Q, P in zip(self.coeff, other.coeff)],
            occ=[cat(q, p) for q, p in zip(self.occ, other.occ)],
            energy=[cat(q, p) for q, p in zip(self.energy, other.energy)],
            idx=[cat(q, p) for q, p in zip(self.idx, other.idx)])

    def __radd__(self, other):
        return self + other

    def __mul__(self, other):
        if not isinstance(other, int):
            raise TypeError("Can't multiply Orbfile by non-int.")
        if other < 0:
            raise ValueError("Other has to be >=0")

        def block_times(n, A):
            return block_diag(*(A for _ in range(n)))

        def cat_times(n, v):
            return np.concatenate([v for _ in range(n)], axis=0)

        n = other
        return self.__class__(
            orbs=[n * n_orbs for n_orbs in self.orbs],
            coeff=[block_times(n, Q) for Q in self.coeff],
            occ=[cat_times(n, q) for q in self.occ],
            energy=[cat_times(n, q) for q in self.energy],
            idx=[cat_times(n, q) for q in self.idx])

    def __rmul__(self, other):
        return self * other

    def __eq__(self, other):
        if not self._same_dims(other):
            return False
        return all([
            (Q == P for Q, P in zip(self.coeff, other.coeff)),
            (q == p for q, p in zip(self.occ, other.occ)),
            (q == p for q, p in zip(self.energy, other.energy)),
            (q == p for q, p in zip(self.idx, other.idx))])

    def __neq__(self, other):
        return not self == other

    def _same_dims(self, other):
        if len(self.orbs) != len(other.orbs):
            return False
        if any(n_s != n_o for n_s, n_o in zip(self.orbs, other.orbs)):
            return False
        return True


    def assimilate(self, other, blocks):
        """Make other blockwise similar to self

    Parameters
    ----------

    other : SpatialOrb
    blocks : list of list of indices
        `blocks[i][j]` contains the index of the j-th block in the i-th irrep
        """
        if not self._same_dims(other):
            raise ValueError("Self and Other have to be same symmetry group "
                             "with the same number of orbitals per irrep")
        new_other = other.copy()
        # Defines a shorter function name
        f = get_kabsch_transformation
        for irrep, indices in enumerate(blocks):
            for GAS_space in indices:
                R = f(self.coeff[irrep][:, GAS_space],
                      other.coeff[irrep][:, GAS_space])
                new_other.coeff[irrep][:, GAS_space] = (
                    other.coeff[irrep][:, GAS_space] @ R.T)
        return new_other


def cat(*args):
    return np.concatenate(args, axis=0)
