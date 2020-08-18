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

from typing import Sequence, Tuple, TypeVar, Type, TextIO
import re
from os import PathLike
from enum import Enum
from copy import deepcopy
from abc import ABCMeta, abstractmethod

import numpy as np
from numpy import array


class FileFormat(Exception):
    pass


class RasOrb_version(Enum):
    v2_0, v2_2 = range(2)


T = TypeVar('T', bound='_Orbitals')


class _Orbitals(metaclass=ABCMeta):

    def __init__(self, orbs, coeff, occ, energy, idx):
        self.orbs = orbs
        self.coeff = coeff
        self.occ = occ
        self.energy = energy
        self.idx = idx

    @classmethod
    @abstractmethod
    def read_orbfile(cls: Type[T], path: PathLike) -> T:
        pass

    @abstractmethod
    def reindex(self,
                new_idx: Sequence[Sequence[int]], inplace: bool=False) -> T:
        pass

    def get_index(self) -> Sequence[Sequence[int]]:
        return [np.arange(n_orbs) for n_orbs in self.orbs]

    def copy(self) -> T:
        return self.__class__(
            orbs=deepcopy(self.orbs),
            coeff=deepcopy(self.coeff),
            occ=deepcopy(self.occ),
            energy=deepcopy(self.energy),
            idx=deepcopy(self.idx))

    def _write_orbfile_stream(self):
        yield from self._write_header()
        yield from self._write_MO_coeff()
        yield from self._write_MO_occ()
        yield from self._write_MO_human_occ()
        yield from self._write_MO_E()
        yield from self._write_CAS_idx()

    def write_orbfile(self, path=None):
        "Write an orbital file. By default an iterator of lines is returned."
        if path is None:
            return self._write_orbfile_stream()
        else:
            with open(path, 'w') as f:
                for line in self._write_orbfile_stream():
                    print(line, file=f)

    @abstractmethod
    def _write_header(self):
        pass

    @abstractmethod
    def _write_MO_coeff(self):
        pass

    @abstractmethod
    def _write_MO_occ(self):
        pass

    @abstractmethod
    def _write_MO_human_occ(self):
        pass

    @abstractmethod
    def _write_MO_E(self):
        pass

    def _write_CAS_idx(self):
        def res_out(v):
            lines = _reshape_output(v, 10, '1')
            return (f'{i % 10} {a}'for i, a in enumerate(lines))

        lines = ['#INDEX']
        for irrep, n_orbs in enumerate(self.orbs):
            lines.append('* 1234567890')
            lines.extend(res_out(self.idx[irrep]))
        return lines

    def __repr__(self):
        return ("{} with {} orbitals per irrep"
                .format(type(self).__name__, self.orbs))


def _reindex(indices: Sequence[Sequence],
             new_indices: Sequence[Sequence[int]]) -> Sequence[Sequence]:
    return [idx[new_idx] for new_idx, idx in zip(new_indices, indices)]


def _spin_reindex(values: Sequence[Sequence],
                  new_indices: Sequence[Sequence[int]]) -> Sequence[Sequence]:
    return {
        spin:
            [v[new_idx] for new_idx, v in zip(new_indices, spin_values)]
        for spin, spin_values in values.items()}


def _read_header(f: TextIO) -> Tuple[RasOrb_version, Sequence[int], bool]:
    line = f.readline()
    while line:
        if '#INPORB' in line:
            version = _read_version(line)
        elif '#INFO' in line:
            spin_orbs, orbs = _read_info(f)
            break
        line = f.readline()
    return version, orbs, spin_orbs


def _read_spin_orbs(
        f: TextIO, orbs: Sequence[int],
        version: RasOrb_version) -> Tuple[Sequence[array], ...]:
    MO_coeff, MO_occ, MO_E = {}, {}, {}
    line = f.readline()
    while line:
        if '#ORB' in line:
            MO_coeff['a'] = _read_coeffs(f, orbs, 5)
        elif '#UORB' in line:
            MO_coeff['b'] = _read_coeffs(f, orbs, 5)
        elif '#OCC' in line:
            MO_occ['a'] = _read_orb_property(
                f, orbs, 10 if version == RasOrb_version.v2_0 else 5)
        elif '#UOCC' in line:
            MO_occ['b'] = _read_orb_property(
                f, orbs, 10 if version == RasOrb_version.v2_0 else 5)
        elif '#ONE' in line:
            MO_E['a'] = _read_orb_property(f, orbs, 10)
        elif '#UONE' in line:
            MO_E['b'] = _read_orb_property(f, orbs, 10)
        elif '#INDEX' in line:
            idx = _read_CAS_idx(f, orbs, 10)
        line = f.readline()
    return MO_coeff, MO_occ, MO_E, idx


def _read_spatial_orbs(
        f: TextIO, orbs: Sequence[int],
        version: RasOrb_version) -> Tuple[Sequence[array], ...]:
    line = f.readline()
    while line:
        if '#ORB' in line:
            MO_coeff = _read_coeffs(f, orbs, 5)
        elif '#OCC' in line:
            MO_occ = _read_orb_property(
                f, orbs, 10 if version == RasOrb_version.v2_0 else 5)
        elif '#ONE' in line:
            MO_E = _read_orb_property(f, orbs, 10)
        elif '#INDEX' in line:
            idx = _read_CAS_idx(f, orbs, 10)
        line = f.readline()
    return MO_coeff, MO_occ, MO_E, idx


def _read_version(line: str) -> RasOrb_version:
    match = re.match('#INPORB (2).([02])', line)
    if match:
        major, minor = match.group(1, 2)
        if major == '2':
            if minor == '0':
                return RasOrb_version.v2_0
            if minor == '2':
                return RasOrb_version.v2_2
        raise ValueError(f'Unknown version {major}.{minor}')
    else:
        raise ValueError('Inporb version could not be determined.')


def _read_info(f: TextIO) -> Tuple[bool, Sequence[int]]:
    _forward(f, 1)
    line = f.readline()
    spin_orbs, n_sym, _ = (int(x) for x in line.split())
    line = f.readline()
    return bool(spin_orbs), [int(irrep) for irrep in line.split()[:n_sym]]


def _read_coeffs(
        f: TextIO, orbs: Sequence[int], cols: int) -> Sequence[array]:
    MO_coeff = [np.empty((n_orbs, n_orbs)) for n_orbs in orbs]
    for irrep, n_orbs in enumerate(orbs):
        rows = (n_orbs + cols - 1) // cols
        for orb in range(n_orbs):
            values = []
            _forward(f, 1)
            for _ in range(rows):
                line = f.readline()
                values.extend([float(x) for x in line.split()])
            MO_coeff[irrep][:, orb] = values
    return MO_coeff


def _read_orb_property(
        f: TextIO, orbs: Sequence[int], cols: int) -> Sequence[array]:
    MO_v = [np.empty(n_orbs) for n_orbs in orbs]
    f.readline()
    for irrep, n_orbs in enumerate(orbs):
        rows = (n_orbs + cols - 1) // cols
        values = []
        for _ in range(rows):
            line = f.readline()
            values.extend([float(x) for x in line.split()])
        MO_v[irrep][:] = values
    return MO_v


def _read_CAS_idx(
        f: TextIO, orbs: Sequence[int], cols: int) -> Sequence[array]:
    idx = []
    for n_orbs in orbs:
        rows = (n_orbs + cols - 1) // cols
        values = []
        next(f)
        for _ in range(rows):
            line = f.readline()
            values.extend(line.strip()[2:])
        idx.append(array(values))
    return idx


def _reshape_output(v, cols, formatter):
    rows = (len(v) + cols - 1) // cols
    for i in range(rows - 1):
        current = v[(i * cols):((i + 1) * cols)]
        yield ''.join(f'{x:{formatter}}' for x in current)
    current = v[(rows - 1) * cols:]
    yield ''.join(f'{x:{formatter}}' for x in current)


def _write_section_1D(title, orbs, values,
                      cols=5, subtitle=None, fmt='22.14E'):
    yield title
    if subtitle is not None:
        yield subtitle
    for irrep in range(len(orbs)):
        yield from _reshape_output(values[irrep][:], cols, fmt)


def _write_section_2D(title, orbs, values, cols=5, subtitle=None,
                      paragraph=None, fmt='22.14E'):
    yield title
    if subtitle is not None:
        yield subtitle
    for irrep, n_orbs in enumerate(orbs):
        for orb in range(n_orbs):
            if paragraph:
                yield paragraph(irrep, orb)
            yield from _reshape_output(values[irrep][:, orb], cols, fmt)


def _forward(f: TextIO, n: int) -> str:
    for _ in range(n):
        line = f.readline()
    return line
