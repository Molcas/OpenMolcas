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

from src.as_methods.base import Table, filter_to_ground
from src.calc_base import SingleRefCalc, MultiRefCalc
from src.utils import GDM_AS


class GDMTable(Table):
    """
    Builds a table displaying absolute ground-state dipole differences by scalar magnitude, where:
    rows = candidate active spaces, columns = reference sources, cells = absolute dipole magnitude difference of ground states.

    :param mr_calcs: list of candidate multiconfigurational active space calculations (of type MultiRefCalc).
    :param reference_sources: mapping of source name (ie, PBE) to a list of SingleRefCalc (derived from reference source input files).
                              in this case (ground state), there should be a single SingleRefCalc per source to represent the ground state calculation.
    """
    PROTOCOL_NAME = GDM_AS

    def __init__(self, mr_calcs: List[MultiRefCalc],
                reference_sources: Dict[str, List[SingleRefCalc]],):
        super().__init__(mr_calcs, reference_sources)


    def calc_diff(self, mr_calc: MultiRefCalc, ref_calcs: List[SingleRefCalc]) -> float:
        """
        Returns the absolute ground-state dipole scalar magnitude difference.
        """

        mr_grounds, ref_grounds = filter_to_ground(mr_calc, ref_calcs)

        if len(ref_grounds) != 1:
            raise ValueError(f'expected exactly one ground-state reference, got {len(ref_grounds)}')
        if len(mr_grounds) != 1:
            raise ValueError(
                f'expected exactly one ground-state entry, got {len(mr_grounds)} in {mr_calc.path}')

        # absolute difference of the multiconfigurational ground state and the reference ground state
        return abs(mr_grounds[0].dipole - ref_grounds[0].state.dipole)
