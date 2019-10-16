from __future__ import annotations
from typing import Sequence, Tuple, TypeVar, Type, TextIO
import re
from os import PathLike
from copy import deepcopy

from numpy import argsort, isclose
from attr import attrs, attrib

from analyze_molcas._base_orbfile import (
    _Orbitals, FileFormat, _reindex,
    _read_spatial_orbs, _read_header, _reshape_output,
    _write_section_1D, _write_section_2D)

from analyze_molcas.spin_orbfile import SpinOrbs


class SpatialOrbs(_Orbitals):

    @classmethod
    def read_orbfile(cls, path: PathLike) -> SpatialOrbs:
        with open(path, 'r') as f:
            version, orbs, spin_orbs = _read_header(f)

            if spin_orbs:
                raise FileFormat('File represents spin orbitals')

            coeff, occ, energy, idx = _read_spatial_orbs(f, orbs, version)
        return cls(orbs, coeff, occ, energy, idx)

    def reindex(self,
            new_idx: Sequence[Sequence[int]], inplace: bool=False) -> SpatialOrbs:
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

    def round_occ(self, inplace: bool=False) -> SpatialOrbs:
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
