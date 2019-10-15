from __future__ import annotations
from typing import Sequence, Tuple, TypeVar, Type, TextIO
import re
from os import PathLike
from enum import Enum
from abc import ABCMeta, abstractmethod
from copy import deepcopy

import numpy as np
from numpy import array, argsort, isclose
from attr import attrs, attrib


class FileFormat(Exception):
    pass


class RasOrb_version(Enum):
    v2_0, v2_2 = range(2)


T = TypeVar('T', bound='_Orbitals')


@attrs
class _Orbitals(metaclass=ABCMeta):
    orbs = attrib()
    coeff = attrib()
    occ = attrib()
    energy = attrib()
    idx = attrib()

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

    def write_orbfile(self, path):
        with open(path, 'w') as f:
            for line in self._get_header():
                print(line, file=f)
            for line in self._write_MO_coeff():
                print(line, file=f)
            for line in self._write_MO_occ():
                print(line, file=f)
            for line in self._write_MO_E():
                print(line, file=f)
            for line in self._write_CAS_idx():
                print(line, file=f)

    def copy(self) -> T:
        return self.__class__(
            orbs=self.orbs.copy(),
            coeff=self.coeff.copy(),
            occ=self.occ.copy(),
            energy=self.energy.copy(),
            idx=self.idx.copy())


class RasOrb(_Orbitals):

    @classmethod
    def read_orbfile(cls, path: PathLike) -> RasOrb:
        with open(path, 'r') as f:
            version, orbs, spin_orbs = _read_header(f)

            if spin_orbs:
                raise FileFormat('File represents spin orbitals')

            coeff, occ, energy, idx = _read_spatial_orbs(f, orbs, version)
        return cls(orbs, coeff, occ, energy, idx)

    def reindex(self,
            new_idx: Sequence[Sequence[int]], inplace: bool=False) -> RasOrb:
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

    def round_occ(self, inplace: bool=False) -> RasOrb:
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

    def to_UhfOrb(self) -> UhfOrb:
        spat = self.round_occ()
        spat.reindex([argsort(-v, kind='stable') for v in spat.occ], True)

        new_occ = {'a': [occ.clip(0, 1) for occ in spat.occ],
                   'b': [(occ - 1).clip(0, 1) for occ in spat.occ]}

        return UhfOrb(
            orbs=spat.orbs.copy(), coeff=spat._spin_copy(spat.coeff),
            occ=new_occ, energy=spat._spin_copy(spat.energy),
            idx=spat.idx.copy())

    @staticmethod
    def _spin_copy(v):
        return {'a': deepcopy(v), 'b': deepcopy(v)}


class UhfOrb(_Orbitals):

    @classmethod
    def read_orbfile(cls, path: PathLike) -> UhfOrb:
        with open(path, 'r') as f:
            version, orbs, spin_orbs = _read_header(f)

            if not spin_orbs:
                raise FileFormat('File represents spatial orbitals')

            coeff, occ, energy, idx = _read_spin_orbs(f, orbs, version)
        return cls(orbs, coeff, occ, energy, idx)

    def reindex(self,
            new_idx: Sequence[Sequence[int]], inplace: bool=False) -> UhfOrb:
        if inplace:
            self.coeff = {
                spin:
                    [coeff[:, idx] for idx, coeff in zip(new_idx, spin_values)]
                for spin, spin_values in self.coeff.items()}
            self.energy = _spin_reindex(self.energy, new_idx)
            self.idx = _reindex(self.idx, new_idx)
            self.occ = _spin_reindex(self.occ, new_idx)
        else:
            new = self.copy()
            new.reindex(new_idx, inplace=True)
            return new



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


def _forward(f: TextIO, n: int) -> str:
    for _ in range(n):
        line = f.readline()
    return line
