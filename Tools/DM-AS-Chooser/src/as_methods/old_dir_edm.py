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

from src.as_methods.base import Table, filter_to_excited
from src.calc_base import SingleRefCalc, MultiRefCalc
from src.utils import OLD_DIR_EDM_AS
import numpy as np


class OldDirEDMTable(Table):
    """
    Builds a table displaying excited-state dipole vector difference:
    rows = candidate active spaces, columns = reference sources, cells = dipole difference of the excited states.

    :param mr_calcs: list of candidate multiconfigurational active space calculations (of type MultiRefCalc).
    :param reference_sources: mapping of source name (ie, PBE) to a list of SingleRefCalc (derived from reference source input files).
    """
    PROTOCOL_NAME = OLD_DIR_EDM_AS

    def __init__(self, mr_calcs: List[MultiRefCalc],
                reference_sources: Dict[str, List[SingleRefCalc]]):
        super().__init__(mr_calcs, reference_sources)

    def calc_diff(self, mr_calc: MultiRefCalc,
                  ref_calcs: List[SingleRefCalc]) -> float:
        """
        Returns the norm of the excited-state dipole vector difference between the candidate and the reference states.
        """
        mr_excited, ref_excited = filter_to_excited(mr_calc, ref_calcs)

        mr_by_state = {s.state: s.vector for s in mr_excited if s.vector is not None}
        ref_by_state = {r.state.state: r.state.vector for r in ref_excited if r.state.vector is not None}

        missing = set(mr_by_state) ^ set(ref_by_state)
        if missing:
            raise ValueError(f'state mismatch between candidate and reference: {missing}')

        if not mr_by_state:
            raise ValueError('no excited-state dipole vectors found among candidate states')

        states = sorted(mr_by_state)
        mr_matrix = np.array([[mr_by_state[s].x, mr_by_state[s].y, mr_by_state[s].z] for s in states])
        ref_matrix = np.array([[ref_by_state[s].x, ref_by_state[s].y, ref_by_state[s].z] for s in states])

        #return float(np.max(np.linalg.norm(mr_matrix - ref_matrix, axis=1)))
        diff_matrix = np.empty(len(mr_matrix))
        for state in range(len(mr_matrix)):
            diff_matrix[state]=abs(np.dot(mr_matrix[state],ref_matrix[state])/np.linalg.norm(ref_matrix[state])-np.linalg.norm(ref_matrix[state]))
        return float(np.max(diff_matrix))
