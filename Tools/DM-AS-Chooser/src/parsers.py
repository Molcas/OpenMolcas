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

"""
Delegates to the molextract package to parse OpenMolcas and ORCA files.

molextract is not a requirement if the user is only providing .csv files.
"""

from typing import List, Dict, Tuple, Optional
import csv

from src.utils import NO_MOLEXTRACT_ERR, DipoleNotFoundError, ElectronicState, NUM_ACTIVE_ELECTRONS_RE, \
    NUM_ACTIVE_ORBITALS_RE, Vector, logger

_RASSCF_MOL_PROPS_RULE = None


def parse_active_space(text: str) -> Tuple[Optional[int], Optional[int]]:
    """
    Reads the active space size directly from the OpenMolcas RASSCF log's "Wave function specifications" summary.

    :param text: the OpenMolcas log.

    :returns: tuple of (# active electrons, # active orbitals). Either value can be None if not found.
    """
    electrons_match = NUM_ACTIVE_ELECTRONS_RE.search(text)
    orbitals_match = NUM_ACTIVE_ORBITALS_RE.search(text)

    num_electrons = int(electrons_match.group(1)) if electrons_match else None
    num_orbitals = int(orbitals_match.group(1)) if orbitals_match else None

    return num_electrons, num_orbitals


def get_orca_parser():
    """
    Sources the ORCADipoleMoment class from molextract package to parse an ORCA .out/.log file.
    """
    try:
        from molextract.parser import Parser
        from molextract.rules.orca.general import ORCADipoleMoment
    except ModuleNotFoundError:
        raise ValueError(NO_MOLEXTRACT_ERR)

    return Parser(ORCADipoleMoment())


def parse_orca_log(path: str) -> List[ElectronicState]:
    """
    Parses the dipole moment block(s) from an ORCA .out/.log file.

    A ground-state-only calculation has exactly one block (SCF, state=0).
    An excited-state calculation has two: SCF first (incidental, state=0), then excited-state second (the requested state).
    """
    parser = get_orca_parser()
    with open(path, 'r') as f:
        raw_log = f.read()

    try:
        data = parser.feed(raw_log)
    except ValueError:
        data = None

    if data is None or len(data) == 0:
        raise DipoleNotFoundError(f'no dipole(s) for {path}')

    return [ElectronicState(state=d['state'],
                            dipole=d['total'],
                            vector=Vector.from_values(d.get('x'), d.get('y'), d.get('z'))) for d in data]


def parse_csv_states(path: str) -> List[ElectronicState]:
    """
     Reads all rows from a CSV with explicit 'state', 'dipole', and 'vector' columns in file order,
      e.g.:
        state, dipole, x, y, z
        0, d0, x0, y0, z0
        1, d1, x1, y1, z1
    """
    with open(path, 'r') as f:
        reader = csv.DictReader(f, skipinitialspace=True)
        rows = list(reader)

    if not rows:
        raise DipoleNotFoundError(f'no dipole(s) for {path}')

    states = []
    for row in rows:
        dipole = None if row.get('dipole') in (None, '') else float(row.get('dipole'))
        vector = Vector.from_values(row.get('x'), row.get('y'), row.get('z'))

        if dipole is None and vector is None:
            logger.warning(f"{path}, state {row['state']}: no dipole or vector provided, skipping")
            continue

        states.append(ElectronicState(state=int(row['state']), dipole=dipole, vector=vector))

    return states


def get_mr_parser():
    """
    Lazily import molextract and return a new parser to parse RASSCF Mol Props
    """
    global _RASSCF_MOL_PROPS_RULE
    try:
        from molextract.rules.molcas import log, general
        from molextract.parser import Parser
    except ModuleNotFoundError:
        raise ValueError(NO_MOLEXTRACT_ERR)

    if _RASSCF_MOL_PROPS_RULE is None:

        class RASSCFMolProps(log.ModuleRule):
            def __init__(self):
                super().__init__('rasscf', [general.MolProps()])

            def reset(self):
                return self.rules[0].reset()

        _RASSCF_MOL_PROPS_RULE = RASSCFMolProps

    return Parser(_RASSCF_MOL_PROPS_RULE())



def parse_mr_log(path: str) -> List[Dict[str, float]]:
    """
    Sources the RASSCFMolProps class from molextract package to parse an OpenMolcas .log file.
    """
    parser = get_mr_parser()
    with open(path, 'r') as f:
        raw_log = f.read()

    try:
        data = parser.feed(raw_log)
    except ValueError:
        # Issue when attempting to parse
        data = None

    if data is None or len(data) == 0:
        raise DipoleNotFoundError(f'no dipole(s) for {path}')

    return data