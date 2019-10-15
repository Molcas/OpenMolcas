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


from analyze_molcas._base_orbfile import (
    _Orbitals, FileFormat, _reindex, _spin_reindex,
    _read_spin_orbs, _read_header)


class SpinOrbs(_Orbitals):

    @classmethod
    def read_orbfile(cls, path: PathLike) -> SpinOrbs:
        with open(path, 'r') as f:
            version, orbs, spin_orbs = _read_header(f)

            if not spin_orbs:
                raise FileFormat('File represents spatial orbitals.')

            coeff, occ, energy, idx = _read_spin_orbs(f, orbs, version)
        return cls(orbs, coeff, occ, energy, idx)

    def reindex(self,
            new_idx: Sequence[Sequence[int]], inplace: bool=False) -> SpinOrbs:
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
