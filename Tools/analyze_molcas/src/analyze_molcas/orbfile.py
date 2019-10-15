from __future__ import annotations
from typing import Sequence, Tuple
import re
from enum import Enum
import numpy as np
from numpy.linalg import svd

from abc import ABCMeta, abstractmethod

from attr import attrs, attrib


class FileFormat(Exception):
    pass


class RasOrb_version(Enum):
    v2_0, v2_2 = range(2)


@attrs
class Orbitals(metaclass=ABCMeta):
    orbs = attrib()
    coeff = attrib()
    occ = attrib()
    energy = attrib()
    idx = attrib()

    @classmethod
    @abstractmethod
    def read_orbfile(cls, path):
        pass

    @abstractmethod
    def reindex(self,
            new_idx: Sequence[Sequence[int]], inplace: bool=False) -> RasOrb:
        pass

    def get_index(self):
        return [np.arange(len(energy)) for energy in self.energy]

    def copy(self):
        return self.__class__(
            orbs=self.orbs.copy(),
            coeff=self.coeff.copy(),
            occ=self.occ.copy(),
            energy=self.energy.copy(),
            idx=self.idx.copy())


class RasOrb(Orbitals):

    @classmethod
    def read_orbfile(cls, path):
        with open(path, 'r') as f:
            version, orbs, spin_orbs = _read_header(f)

            if spin_orbs:
                raise FileFormat('File represents spin orbitals')

            coeff, occ, energy, idx = _read_spin_orbs(f, orbs, version)
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


def _reindex(indices: Sequence[Sequence],
             new_indices: Sequence[Sequence[int]]) -> Sequence[Sequence]:
    return [idx[new_idx] for new_idx, idx in zip(new_indices, indices)]



class UhfOrb(Orbitals):

    @classmethod
    def read_orbfile(cls, path):
        with open(path, 'r') as f:
            version, orbs, spin_orbs = _read_header(f)

            if not spin_orbs:
                raise FileFormat('File represents spatial orbitals')

            coeff, occ, energy, idx = _read_spin_orbs(f, orbs, version)
        return cls(orbs, coeff, occ, energy, idx)


    def reindex(self,
            new_idx: Sequence[Sequence[int]], inplace: bool=False) -> RasOrb:
        if inplace:
            self.coeff = {
                spin: [coeff[:, idx] for idx, coeff in zip(new_idx, coeffs)]
                for spin, coeffs in self.coeff.items()}
            self.energy = [
                energy[idx] for idx, energy in zip(new_idx, self.energy)]
            self.idx = [kind[idx] for idx, kind in zip(new_idx, self.idx)]
            self.occ = [occ[idx] for idx, occ in zip(new_idx, self.occ)]
            self.orbs = [len(idx) for idx in new_idx]
        else:
            new = self.copy()
            new.reindex(new_idx, inplace=True)
            return new


def _read_header(f):
    line = f.readline()
    while line:
        if '#INPORB' in line:
            version = _read_version(line)
        elif '#INFO' in line:
            spin_orbs, orbs = _read_info(f)
            break
        line = f.readline()
    return version, orbs, spin_orbs


def _read_spin_orbs(f, orbs, version):
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


def _read_spatial_orbs(f, orbs, version):
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


def _read_info(f) -> Tuple[bool, Sequence[int]]:
    _forward(f, 1)
    line = f.readline()
    spin_orbs, n_sym, _ = (int(x) for x in line.split())
    line = f.readline()
    return bool(spin_orbs), [int(irrep) for irrep in line.split()[:n_sym]]


def _read_coeffs(f, orbs, cols):
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


def _read_orb_property(f, orbs, cols):
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


def _read_CAS_idx(f, orbs, cols):
    idx = []
    for n_orbs in orbs:
        rows = (n_orbs + cols - 1) // cols
        values = []
        next(f)
        for _ in range(rows):
            line = f.readline()
            values.extend(line.strip()[2:])
        idx.append(np.array(values))
    return idx


def _forward(f, n):
    for _ in range(n):
        line = f.readline()
    return line
