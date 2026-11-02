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
from src.utils import EDM_AS
import numpy as np


class EDMTable(Table):
    """
    Builds a table displaying max absolute excited-state dipole differences by scalar magnitude, where:
    rows = candidate active spaces, columns = reference sources, cells = max absolute dipole magnitude difference across all excited states.

    :param mr_calcs: list of candidate multiconfigurational active space calculations (of type MultiRefCalc).
    :param reference_sources: mapping of source name (ie, PBE) to a list of SingleRefCalc (derived from reference source input files).
    """
    PROTOCOL_NAME = EDM_AS

    def __init__(self, mr_calcs: List[MultiRefCalc],
                reference_sources: Dict[str, List[SingleRefCalc]]):
        super().__init__(mr_calcs, reference_sources)


    def calc_diff(self, mr_calc: MultiRefCalc,
                      ref_calcs: List[SingleRefCalc]) -> float:
            """
            Returns the max absolute excited-state dipole scalar magnitude difference.
            """
            mr_excited, ref_excited = filter_to_excited(mr_calc, ref_calcs)

            # {excited_state_val: magnitude}
            mr_by_state = {s.state: s.dipole for s in mr_excited if s.dipole is not None}
            ref_by_state = {r.state.state: r.state.dipole for r in ref_excited if r.state.dipole is not None}

            missing = set(mr_by_state) ^ set(ref_by_state)
            if missing:
                raise ValueError(f'state mismatch between candidate and reference: {missing}')

            if not mr_by_state:
                raise ValueError('no excited-state dipole values found among candidate states')

            # Sorted array of magnitudes
            states = sorted(mr_by_state)
            mr_arr = np.array([mr_by_state[s] for s in states])
            ref_arr = np.array([ref_by_state[s] for s in states])

            # Maximum absolute difference of multiconfiguration and reference states.
            # State order ensured by previous sorting
            return float(np.max(np.abs(mr_arr - ref_arr)))