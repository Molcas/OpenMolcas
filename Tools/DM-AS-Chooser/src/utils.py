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


import re
import os
import logging
from dataclasses import dataclass
from typing import Optional

# Protocols
GDM_AS = 'gdm-as'
EDM_AS = 'edm-as'
D2DM_AS = 'd2dm-as'
DIR_GDM_AS = 'dir-gdm-as'
DIR_EDM_AS = 'dir-edm-as'
DIR_D2DM_AS = 'dir-d2dm-as'
OLD_DIR_GDM_AS = 'old-dir-gdm-as'
OLD_DIR_EDM_AS = 'old-dir-edm-as'

# Regexes
ACTIVE_SPACE_RE = re.compile(r'(\d+)-(\d+)')
NUM_ACTIVE_ELECTRONS_RE = re.compile(r'Number of electrons in active shells\s+(\d+)')
NUM_ACTIVE_ORBITALS_RE = re.compile(r'Number of active orbitals\s+(\d+)')

# File extensions
CSV_EXT = '.csv'
LOG_EXT = '.log'
ORCA_EXT = '.orca'
OUT_EXT = '.out'
VALID_EXTS = (CSV_EXT, LOG_EXT, ORCA_EXT, OUT_EXT)
VALID_ORCA_EXTS = (ORCA_EXT, OUT_EXT, LOG_EXT)

# Error messages
NO_MOLEXTRACT_ERR = 'you must have molextract installed to parse log files'

# Logging
logger = logging.getLogger('active_space_chooser')
handler = logging.StreamHandler()
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)

# Conversion Units
CONSTANTS = {
    'DEBYE_PER_AU': 2.541765,
    'HARTREE_TO_EV': 27.211386245981,
}

# Ground state discrepancy tolerance
GROUND_STATE_TOL = 1e-4

# Calculations the extraction scripts can be run on
EXTRACTABLE_KEYWORDS = {'CASSCF', 'CASPT2', 'CAS-PDFT'}


class DipoleNotFoundError(Exception):
    """
    To represent an error that occurs if a parser cannot locate the dipole moment in a given file.
    """
    pass

@dataclass
class Vector:
    x: float
    y: float
    z: float

    @classmethod
    def from_values(cls, x, y, z) -> Optional['Vector']:
        """
        Builds a Vector from three raw component values.
        """
        if x is None or y is None or z is None or x == '' or y == '' or z == '':
            return None
        return cls(x=float(x), y=float(y), z=float(z))

    def __eq__(self, other: 'Vector') -> bool:
        if not isinstance(other, Vector):
            return False
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __hash__(self) -> int:
        return hash((self.x, self.y, self.z))

@dataclass
class ElectronicState:
    """
    An electronic state and its dipole moment.

    :param state: a state of 0 represents the ground state and 1+ represents an excited state.
    :param dipole: the dipole moment of the electronic state.
    :param vector: the dipole moment vector
    """
    state: int
    dipole: Optional[float] = None
    vector: Optional[Vector] = None

    def __str__(self) -> str:
        """
        Returns a human-readable summary of this state.
        """
        lines = [f'state = {self.state}'
                 f'dipole = {self.dipole}'
                 f'vector = {self.vector}']

        return '\n'.join(lines)

    def __eq__(self, other: 'ElectronicState') -> bool:
        if not isinstance(other, ElectronicState):
            return False
        return self.state == other.state and self.dipole == other.dipole and self.vector == other.vector

    def __hash__(self) -> int:
        return hash((self.state, self.dipole, self.vector))

def check_valid_file(path: str):
    """
    To validate that the given file path exists and has an accepted file extension.
    """
    if not os.path.exists(path):
        raise ValueError(f'{path} does not exist')
    _, ext = os.path.splitext(path)
    if ext not in VALID_EXTS:
        raise ValueError(f'{path} must be one of {VALID_EXTS}')







