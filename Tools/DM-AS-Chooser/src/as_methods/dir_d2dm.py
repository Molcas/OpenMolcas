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
from src.utils import DIR_D2DM_AS
import numpy as np


class DirD2DMTable(Table):
    """
    Builds a table displaying absolute dipole differences by vectors of states of interest where:
    rows = candidate active spaces, columns = reference sources, cells = dipole difference of the states of interest.

    States of interest are the states with the largest absolute difference between ground and excited states.

    :param mr_calcs: list of candidate multiconfigurational active space calculations (of type MultiRefCalc).
    :param reference_sources: mapping of source name (ie, PBE) to a list of SingleRefCalc (derived from reference source input files).
    """
    PROTOCOL_NAME = DIR_D2DM_AS

    def __init__(self, mr_calcs: List[MultiRefCalc],
                reference_sources: Dict[str, List[SingleRefCalc]]):
        super().__init__(mr_calcs, reference_sources)

    def calc_diff(self, mr_calc: MultiRefCalc,
                  ref_calcs: List[SingleRefCalc]) -> float:
        """
        Compares dipoles via vector-norm difference (not magnitude).

        Steps:
            1. Get the ground state and excited state vectors.
            2. Create an array of (excited state vectors - ground state vector).
            3. Get the index of the largest (excited state vector - ground state vector) value = state of interest.
                Should be two SOIs: one for reference and one for multiconfigurational candidates.
            4. Final value is the normalized difference of the two SOIs.
        """

        # Ground State References
        mr_ground = next((s for s in mr_calc.states if s.state == 0), None)
        ref_ground = next((s for s in ref_calcs if s.state.state == 0), None)
        if mr_ground is None or mr_ground.vector is None:
            raise ValueError('no ground-state entry found among candidate states')
        if ref_ground is None or ref_ground.state.vector is None:
            raise ValueError('no ground-state entry found among reference states')

        check_consistent_ground_states(ref_calcs)

        # Excited States
        mr_excited = [s for s in mr_calc.states if s.state != 0 and s.vector is not None]
        ref_excited = [s for s in ref_calcs if s.state.state != 0 and s.state.vector is not None]

        # Ground States as Vectors
        mr_ground_vec = np.array([mr_ground.vector.x, mr_ground.vector.y, mr_ground.vector.z])
        ref_ground_vec = np.array([ref_ground.state.vector.x, ref_ground.state.vector.y,
                                   ref_ground.state.vector.z])

        # Excited States as Vectors
        mr_vectors = np.array([[s.vector.x, s.vector.y, s.vector.z] for s in mr_excited])
        ref_vectors = np.array([[s.state.vector.x, s.state.vector.y, s.state.vector.z]
                                for s in ref_excited])

        if mr_vectors.size == 0:
            raise ValueError('no excited-state dipole vectors found among candidate states')
        if ref_vectors.size == 0:
            raise ValueError('no excited-state dipole vectors found among reference states')

        # Excited State Vectors - Ground State Vector
        mr_diffs = np.linalg.norm(mr_vectors - mr_ground_vec, axis=1)
        ref_diffs = np.linalg.norm(ref_vectors - ref_ground_vec, axis=1)

        # Max Vector Difference = SOI
        mr_soi_vec = mr_vectors[np.argmax(mr_diffs)]
        ref_soi_vec = ref_vectors[np.argmax(ref_diffs)]

        # Normalized Difference of SOIs
        return float(np.linalg.norm(mr_soi_vec - ref_soi_vec))