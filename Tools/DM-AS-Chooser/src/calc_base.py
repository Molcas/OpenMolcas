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

from dataclasses import dataclass, field
from typing import List, Optional

from src.parsers import parse_active_space
from src.utils import ElectronicState, check_valid_file, ACTIVE_SPACE_RE, logger, Vector, CONSTANTS


@dataclass
class MultiRefCalc:
    """
    To represent a multiconfigurational calculation.
    These are used by the protocols as the candidate active spaces.

    :param num_electrons: the number of electrons included in the active space.
    :param num_orbitals: the number of molecular orbitals included in the active space.
    :param path: the path to the OpenMolcas .log file or .csv source file.
    :param states: the list of electronic states of the molecule, both ground state and any excited states.
    """
    num_electrons: Optional[int]
    num_orbitals: Optional[int]
    path: str
    states: List[ElectronicState] = field(default_factory=list, init=False)

    def __post_init__(self):
        check_valid_file(self.path)
        self.states = self._parse_states()

    @classmethod
    def infer(cls, path: str) -> 'MultiRefCalc':
        """
        Constructs a MultiRefCalc from a file path.

        If given an OpenMolcas .log file, sources the number of active electrons and orbitals from the file itself.
        Therefore, the name of the OpenMolcas .log file is insignificant.

        If given a .csv file, assumes that the file name contains #electrons-#orbitals.

        :param path: the path to the OpenMolcas .log file or .csv source file.

        :returns: a MultiRefCalc object representing the electronic states of the molecule from the given file path.
        """

        check_valid_file(path)

        if path.endswith('.csv'):
            m = ACTIVE_SPACE_RE.search(path)
            if m:
                num_electrons = int(m.group(1))
                num_orbitals = int(m.group(2))
            else:
                num_electrons = None
                num_orbitals = None
        else:
            with open(path, 'r') as f:
                text = f.read()
            num_electrons, num_orbitals = parse_active_space(text)

        if num_electrons is None or num_orbitals is None:
            logger.warning(f'unable to determine active space for {path}')

        return cls(num_electrons, num_orbitals, path)

    def _parse_states(self) -> List[ElectronicState]:
        """
        Parses the different electronic states of the molecule from the source file.

        :returns: a list of ElectronicState objects representing all electronic states represented in the source file.
        """
        from src.parsers import parse_mr_log, parse_csv_states

        if self.path.endswith('.csv'):
            return parse_csv_states(self.path)
        elif self.path.endswith('.log'):
            data = parse_mr_log(self.path)
            return [ElectronicState(state=i,
                                    dipole=e['dipole']['total'],
                                    vector=Vector(x=e['dipole']['x'],
                                                  y=e['dipole']['y'],
                                                  z=e['dipole']['z']))
                    for i, e in enumerate(data)]
        else:
            raise ValueError(f'unrecognized file extension: {self.path}')

@dataclass
class SingleRefCalc:
    """
    To represent a single reference calculation.
    These are used by the protocols as the reference dipole moments.

    :param state: the electronic state of the molecule.
    :param path: the path to the ORCA .log file or .csv source file.
    """
    state: ElectronicState
    path: str

    @staticmethod
    def _read_states(path: str) -> List[ElectronicState]:
        """
        Reads the given file and converts it's content to a list of ElectronicState objects.

        NOTE: assumes that the orca vector information is in a.u. and automatically converts to Debye.

        :param path: the path to the ORCA .log file or .csv source file.

        :returns: a list of ElectronicState objects representing all electronic states in the source file.
        """
        from src.parsers import parse_orca_log, parse_csv_states

        check_valid_file(path)

        if path.endswith('.csv'):
            return parse_csv_states(path)
        else:
            states = parse_orca_log(path)

            for state in states:
                if state.vector is not None:
                    state.vector = Vector(state.vector.x * CONSTANTS['DEBYE_PER_AU'],
                                          state.vector.y * CONSTANTS['DEBYE_PER_AU'],
                                          state.vector.z * CONSTANTS['DEBYE_PER_AU'])
            return states

    @classmethod
    def all_from_file(cls, path: str) -> List['SingleRefCalc']:
        """
        Reads every state from a reference file (ground and excited alike) as its own SingleRefCalc.

        :param path: the path to the ORCA .log file or .csv source file.

        :returns: a list of SingleRefCalc objects representing all electronic states with their dipoles in the source file.
        """
        states = cls._read_states(path)
        return [cls(state=state, path=path) for state in states]
