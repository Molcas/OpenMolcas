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

from typing import Dict, List

from src.as_methods.base import Table, check_consistent_ground_states
from src.calc_base import SingleRefCalc, MultiRefCalc
from src.utils import D2DM_AS
import numpy as np


class D2DMTable(Table):
    """
    Builds a table displaying absolute dipole differences by scalar magnitude of states of interest where:
    rows = candidate active spaces, columns = reference sources, cells = dipole difference of the states of interest.

    States of interest are the states with the largest absolute difference between ground and excited states.

    :param mr_calcs: list of candidate multiconfigurational active space calculations (of type MultiRefCalc).
    :param reference_sources: mapping of source name (ie, PBE) to a list of SingleRefCalc (derived from reference source input files).
    """
    PROTOCOL_NAME = D2DM_AS

    def __init__(self, mr_calcs: List[MultiRefCalc],
                reference_sources: Dict[str, List[SingleRefCalc]]):
        super().__init__(mr_calcs, reference_sources)

    def calc_diff(self, mr_calc: MultiRefCalc,
                  ref_calcs: List[SingleRefCalc]) -> float:
        """
        Compares dipoles via scalar magnitude differences of states of interest.

        Steps:
            1. Get the ground state and excited state magnitudes.
            2. Create an array of (excited state magnitudes - ground state magnitude).
            3. Get the index of the largest (excited state magnitude - ground state magnitude) value = state of interest.
                Should be two SOIs: one for reference and one for multiconfigurational candidates.
            4. Final value is the absolute difference of the two SOIs.
        """

        # Ground State References
        mr_ground = next((s for s in mr_calc.states if s.state == 0), None)
        ref_ground = next((s for s in ref_calcs if s.state.state == 0), None)
        if mr_ground is None or mr_ground.dipole is None:
            raise ValueError('no ground-state entry found among candidate states')
        if ref_ground is None or ref_ground.state.dipole is None:
            raise ValueError('no ground-state entry found among reference states')

        check_consistent_ground_states(ref_calcs)

        # Excited States
        mr_excited = [s for s in mr_calc.states if s.state != 0]
        ref_excited = [s for s in ref_calcs if s.state.state != 0]

        # Excited State Scalar Magnitudes
        mr_dipoles = np.array([s.dipole for s in mr_excited if s.dipole is not None])
        ref_dipoles = np.array([s.state.dipole for s in ref_excited if s.state.dipole is not None])
        if mr_dipoles.size == 0:
            raise ValueError('no excited-state dipole values found among candidate states')
        if ref_dipoles.size == 0:
            raise ValueError('no excited-state dipole values found among reference states')

        # Excited State Mags - Ground State Mag
        mr_diffs = np.abs(mr_dipoles - mr_ground.dipole)
        ref_diffs = np.abs(ref_dipoles - ref_ground.state.dipole)

        # Max Mag Difference = SOI
        mr_soi_dipole = mr_dipoles[np.argmax(mr_diffs)]
        ref_soi_dipole = ref_dipoles[np.argmax(ref_diffs)]

        # Normalized Difference of SOIs
        return float(abs(mr_soi_dipole - ref_soi_dipole))