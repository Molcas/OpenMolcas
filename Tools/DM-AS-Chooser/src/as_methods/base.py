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
# Copyright (C) 2026, Benjamin W. Kaufold                              *
#               2026, Meredith E. Wanta                                *
#               2026, Sijia S. Dong                                    *
#***********************************************************************


import csv
import math
from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Dict, Tuple
import numpy as np
from src.calc_base import MultiRefCalc, SingleRefCalc
from src.utils import logger, GROUND_STATE_TOL, ElectronicState


def check_consistent_ground_states(ref_calcs: List[SingleRefCalc]) -> None:
    """
    Verifies the ground-state dipole moments agree.

    :param ref_calcs: List of reference calculations to compare.

    :raises ValueError: if any two ground-state entries disagree on magnitude, or (whenever a vector is
        present on both) on direction, beyond GROUND_STATE_TOL.
    """
    grounds = [r for r in ref_calcs if r.state.state == 0]
    if len(grounds) <= 1:
        return
    first = grounds[0]
    for other in grounds[1:]:
        if not math.isclose(first.state.dipole, other.state.dipole, abs_tol=GROUND_STATE_TOL):
            raise ValueError(
                f'inconsistent ground-state dipole among reference files: '
                f'{first.path} gives {first.state.dipole}, {other.path} gives {other.state.dipole}')
        if first.state.vector is not None and other.state.vector is not None:
            first_vec = np.array([first.state.vector.x, first.state.vector.y, first.state.vector.z])
            other_vec = np.array([other.state.vector.x, other.state.vector.y, other.state.vector.z])
            if np.linalg.norm(first_vec - other_vec) > GROUND_STATE_TOL:
                raise ValueError(
                    f'inconsistent ground-state vector among reference files: '
                    f'{first.path} gives {first_vec.tolist()}, {other.path} gives {other_vec.tolist()}')

def filter_to_ground(mr_calc: MultiRefCalc, ref_calcs: List[SingleRefCalc]
                     ) -> Tuple[List[ElectronicState], List[SingleRefCalc]]:
    """
    Filters both sides down to ground-state (state==0) entries only.
    Reference entries are verified for consistency and collapsed to a single representative entry if more than one is found.
    Used by GDMTable and DirGDMTable.

    :returns: (mr_ground_states, ref_ground_calcs), the candidate's ground-state entries (normally
        exactly one) and the reference's, collapsed to at most one.
    """
    mr_ground_states = [s for s in mr_calc.states if s.state == 0]
    ref_ground_calcs = [r for r in ref_calcs if r.state.state == 0]
    check_consistent_ground_states(ref_ground_calcs)
    return mr_ground_states, ref_ground_calcs[:1]

def filter_to_excited(mr_calc: MultiRefCalc, ref_calcs: List[SingleRefCalc]
                      ) -> Tuple[List, List[SingleRefCalc]]:
    """
    Filters both sides down to excited-state (state != 0) entries only.
    Used by EDMTable and DirEDMTable.

    :returns: (mr_excited_states, ref_excited_calcs).
    """
    mr_excited = [s for s in mr_calc.states if s.state != 0]
    ref_excited = [r for r in ref_calcs if r.state.state != 0]
    return mr_excited, ref_excited

class Table(ABC):
    """
    Abstract base class for active space selection methods.

    Each subclass must implement the following:
     - PROTOCOL_NAME: the CLI input string (e.g. "gdm-as") to represent the active space selection method being used.
     - calc_diff: the method used to calculate the difference between the candidate and reference dipole moments.


    :param mr_calcs: (List[MultiRefCalc]) list of candidate multiconfigurational active space calculations.
    :param reference_sources: (Dict[str, List[SingleRefCalc]]) mapping of source name (ie, PBE) to a list of SingleRefCalc
                              (derived from reference source input files).
    """

    PROTOCOL_NAME: str

    def __init__(self, mr_calcs: List[MultiRefCalc],
                reference_sources: Dict[str, List[SingleRefCalc]]):
        self.mr_calcs = mr_calcs
        self.reference_sources = reference_sources

    def __str__(self) -> str:
        """
        Returns a human-readable summary of this table's fields.

        :returns: a multi-line string summarizing self.mr_calcs and self.reference_sources.
        """
        lines = [f'{self.PROTOCOL_NAME}:', f'  {len(self.mr_calcs)} candidate(s):']

        for calc in self.mr_calcs:
            label = f'{calc.num_electrons:02d}-{calc.num_orbitals:02d}'
            lines.append(f'    {label}  ({calc.path})')

        lines.append(f'  {len(self.reference_sources)} reference source(s):')
        for name, ref_calcs in self.reference_sources.items():
            lines.append(f'    {name}: {len(ref_calcs)} file(s)')
            for ref_calc in ref_calcs:
                lines.append(f'      {ref_calc.path}')

        return '\n'.join(lines)


    @abstractmethod
    def calc_diff(self, mr_calc: MultiRefCalc, ref_calcs: List[SingleRefCalc]) -> float:
        """
        Calculates the dipole moment difference relative to the protocol.

        :param mr_calc: MultiRefCalc object, contains all candidate active spaces.
        :param ref_calcs: list of SingleRefCalc objects representing reference sources.

        :returns: float, representing the dipole moment difference relative to the protocol to be written in the output.tsv file.
        """
        raise NotImplementedError

    def build(self) -> Tuple[List[str], List[List]]:
        """
        Builds a table of active space calculations by dipole moment differences (relative to selected protocol).

        :returns: (columns, rows) where columns is the list of reference source names (in column order) and
         rows is a list of [label, val1, val2, ...], one per candidate active space, where val is '' if that cell couldn't be computed.
        """
        columns = list(self.reference_sources.keys())
        rows = []

        for mr_calc in self.mr_calcs:
            label = f'{mr_calc.num_electrons:02d}-{mr_calc.num_orbitals:02d}'
            row = [label]

            for name in columns:
                try:
                    value = self.calc_diff(mr_calc, self.reference_sources[name])
                except ValueError as e:
                    logger.warning(f'{mr_calc.path} / {name}: {e}')
                    value = ''
                row.append(value)

            rows.append(row)

        return columns, rows



    def write_tsv(self, output_path: Path) -> Path:
        """
        Writes the built table as a TSV file.

        :param output_path: path to output TSV file. Defaults to "results/<protocol>.csv".
        :returns: path to output TSV file
        """
        path = output_path / f'{self.PROTOCOL_NAME}.tsv'

        columns, rows = self.build()
        with open(path, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow([''] + columns)
            writer.writerows(rows)

        return path