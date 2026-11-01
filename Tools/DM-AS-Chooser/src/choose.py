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
Picks the winning active space for each reference source.
"""
from typing import Dict, List, Optional, Union


def pick_winners(columns: List[str],
                 rows: List[List[Union[str, float]]]) -> Dict[str, Optional[str]]:
    """
    For each column (reference source), find the row with the smallest
    value and return its label.

    :param columns: column names, as returned by Table.build().
    :param rows: rows as [label, val1, val2, ...], same shape as returned by build().
                 Empty-string cells (failed computations) are skipped.

    :return: Mapping of column name -> winning active-space label.
    """
    winners = {}

    for col_index, column in enumerate(columns, start=1):
        best_label = None
        best_value = None

        for row in rows:
            label = row[0]
            value = row[col_index]

            if value == '':
                continue

            value = float(value)
            if best_value is None or value < best_value:
                best_value = value
                best_label = label

        winners[column] = best_label

    return winners