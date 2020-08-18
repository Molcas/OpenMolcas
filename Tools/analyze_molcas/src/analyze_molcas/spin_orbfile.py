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

from typing import Sequence
from os import PathLike


from analyze_molcas._base_orbfile import (
    _Orbitals, FileFormat, _reindex, _spin_reindex,
    _read_spin_orbs, _read_header,
    _write_section_1D, _write_section_2D)


class SpinOrbs(_Orbitals):

    @classmethod
    def read_orbfile(cls, path: PathLike):
        with open(path, 'r') as f:
            version, orbs, spin_orbs = _read_header(f)

            if not spin_orbs:
                raise FileFormat('File represents spatial orbitals.')

            coeff, occ, energy, idx = _read_spin_orbs(f, orbs, version)
        return cls(orbs, coeff, occ, energy, idx)

    def reindex(self,
                new_idx: Sequence[Sequence[int]], inplace: bool=False):
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

    def _write_header(self):
        return ['#INPORB 2.2',
                '#INFO',
                '* PyOrb written',
                f"{1:8d}{len(self.orbs):8d}{0:8d}",
                f"{''.join(f'{x:8d}' for x in self.orbs)}",
                f"{''.join(f'{x:8d}' for x in self.orbs)}",
                '* SpinOrbs written from analyze_molcas']

    def _write_MO_coeff(self):
        def get_paragraph(irrep, orb):
            return f'* ORBITAL{irrep + 1:5d}{orb + 1:5d}'

        yield from _write_section_2D(
            title='#ORB', paragraph=get_paragraph,
            orbs=self.orbs, values=self.coeff['a'])
        yield from _write_section_2D(
            title='#UORB', paragraph=get_paragraph,
            orbs=self.orbs, values=self.coeff['b'])

    def _write_MO_occ(self):
        yield from _write_section_1D(
            title='#OCC', subtitle='* OCCUPATION NUMBERS',
            orbs=self.orbs, values=self.occ['a'])
        yield from _write_section_1D(
            title='#UOCC', subtitle='* Beta OCCUPATION NUMBERS',
            orbs=self.orbs, values=self.occ['b'])

    def _write_MO_human_occ(self):
        yield from _write_section_1D(
            title='#OCHR', subtitle='* OCCUPATION NUMBERS (HUMAN-READABLE)',
            orbs=self.orbs, values=self.occ['a'], cols=10, fmt='8.4f')
        yield from _write_section_1D(
            title='#UOCHR',
            subtitle='* Beta OCCUPATION NUMBERS (HUMAN-READABLE)',
            orbs=self.orbs, values=self.occ['b'], cols=10, fmt='8.4f')

    def _write_MO_E(self):
        yield from _write_section_1D(
            title='#ONE', subtitle='* ONE ELECTRON ENERGIES',
            orbs=self.orbs, values=self.energy['a'], cols=10, fmt='12.4E')
        yield from _write_section_1D(
            title='#UONE', subtitle='* Beta ONE ELECTRON ENERGIES',
            orbs=self.orbs, values=self.energy['b'], cols=10, fmt='12.4E')
