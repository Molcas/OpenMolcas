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
Tests for choose.py's pick_winners.
"""
from src.choose import pick_winners


class TestPickWinners:
    def test_picks_smallest_value_per_column(self):
        columns = ['BLYP', 'PBE']
        rows = [
            ['06-06', 0.5, 0.2],
            ['06-07', 0.1, 0.9],
            ['08-08', 0.3, 0.05],
        ]

        assert pick_winners(columns, rows) == {'BLYP': '06-07', 'PBE': '08-08'}

    def test_skips_empty_string_cells(self):
        columns = ['BLYP', 'PBE']
        rows = [
            ['06-06', '', 0.2],
            ['06-07', 0.1, ''],
            ['08-08', 0.3, 0.05],
        ]

        assert pick_winners(columns, rows) == {'BLYP': '06-07', 'PBE': '08-08'}

    def test_column_with_no_valid_values_gives_none(self):
        columns = ['BLYP', 'PBE']
        rows = [
            ['06-06', '', 0.2],
            ['06-07', '', 0.9],
        ]

        assert pick_winners(columns, rows) == {'BLYP': None, 'PBE': '06-06'}

    def test_no_rows_at_all_gives_none_for_every_column(self):
        assert pick_winners(['BLYP', 'PBE'], []) == {'BLYP': None, 'PBE': None}

    def test_tie_keeps_the_first_occurrence(self):
        columns = ['BLYP', 'PBE']
        rows = [
            ['06-06', 0.5, 0.1],
            ['06-07', 0.5, 0.1],
        ]

        assert pick_winners(columns, rows) == {'BLYP': '06-06', 'PBE': '06-06'}

    def test_negative_value_beats_smaller_magnitude_positive(self):
        columns = ['BLYP']
        rows = [
            ['06-06', -5.0],
            ['06-07', 0.5],
        ]

        assert pick_winners(columns, rows) == {'BLYP': '06-06'}

    def test_single_row(self):
        assert pick_winners(['BLYP'], [['06-06', 0.5]]) == {'BLYP': '06-06'}

    def test_string_typed_values_are_converted_correctly(self):
        columns = ['BLYP', 'PBE']
        rows = [
            ['06-06', '0.5', '0.2'],
            ['06-07', '0.1', '0.9'],
        ]

        assert pick_winners(columns, rows) == {'BLYP': '06-07', 'PBE': '06-06'}