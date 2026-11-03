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
Tests CSV file parsing.
"""
import pytest

from src.as_methods.dir_gdm import DirGDMTable
from src.as_methods.gdm import GDMTable
from method_tests.mock_utils import FakeMultiRefCalc, FakeSingleRefCalc
from src.parsers import parse_csv_states
from src.utils import ElectronicState, Vector
from unittest.mock import patch, mock_open

CSV_ALL = """state,dipole,x,y,z
0,1.9812,0.0,0.0,-1.9812
1,5.2057,0.15,-0.02,5.1980"""

CSV_DIPOLE = """state,dipole
0,1.9812
1,5.2057"""

CSV_VECTOR = """state,x,y,z
0,0.0,0.0,-1.9812
1,0.15,-0.02,5.1980"""


class TestOptionalFields:
    def test_optional_vector(self):
        result = [ElectronicState(state=0, dipole=1.9812), ElectronicState(state=1, dipole=5.2057)]

        with patch('builtins.open', mock_open(read_data=CSV_DIPOLE)):
            states = parse_csv_states('dipole.csv')

        assert states == result

    def test_optional_dipole(self):
        vec1 = Vector(x=0.0, y=0.0, z=-1.9812)
        vec2 = Vector(x=0.15, y=-0.02, z=5.1980)

        result = [ElectronicState(state=0, vector=vec1), ElectronicState(state=1, vector=vec2)]

        with patch('builtins.open', mock_open(read_data=CSV_VECTOR)):
            states = parse_csv_states('vector.csv')

        assert states == result

    def test_both(self):
        vec1 = Vector(x=0.0, y=0.0, z=-1.9812)
        vec2 = Vector(x=0.15, y=-0.02, z=5.1980)

        result = [ElectronicState(state=0, dipole=1.9812, vector=vec1), ElectronicState(state=1, dipole=5.2057, vector=vec2)]

        with patch('builtins.open', mock_open(read_data=CSV_ALL)):
            states = parse_csv_states('all.csv')

        assert states == result

    def test_optional_vector_runs(self):
        state0 = ElectronicState(state=0, dipole=2.9812)
        state1 = ElectronicState(state=1, dipole=5.0000)

        with patch('builtins.open', mock_open(read_data=CSV_DIPOLE)):
            states = [parse_csv_states('dipole.csv')[0]]

        mr_calc = FakeMultiRefCalc(states=states)
        ref_calcs = [FakeSingleRefCalc(state=state0)]

        table = GDMTable(mr_calcs=[], reference_sources={})
        diff = table.calc_diff(mr_calc, ref_calcs)

        assert diff == pytest.approx(1.000, abs=0.001)

    def test_optional_dipole_runs(self):
        vec1 = Vector(x=0.0, y=0.0, z=-1.9812)
        vec2 = Vector(x=0.15, y=-0.02, z=5.1980)

        state0 = ElectronicState(state=0, vector=vec2)

        with patch('builtins.open', mock_open(read_data=CSV_VECTOR)):
            states = [parse_csv_states('vector.csv')[0]]

        mr_calc = FakeMultiRefCalc(states=states)
        ref_calcs = [FakeSingleRefCalc(state=state0)]

        table = DirGDMTable(mr_calcs=[], reference_sources={})
        diff = table.calc_diff(mr_calc, ref_calcs)

        assert diff == pytest.approx(7.18, abs=0.001)

    def test_both_run(self):
        vec1 = Vector(x=0.0, y=0.0, z=-1.9812)
        state1 = ElectronicState(state=0, dipole=5.0000, vector=vec1)

        with patch('builtins.open', mock_open(read_data=CSV_ALL)):
            states = [parse_csv_states('all.csv')[0]]

        mr_calc = FakeMultiRefCalc(states=states)
        ref_calcs = [FakeSingleRefCalc(state=state1)]

        table_dir = DirGDMTable(mr_calcs=[], reference_sources={})
        diff1 = table_dir.calc_diff(mr_calc, ref_calcs)

        table = GDMTable(mr_calcs=[], reference_sources={})
        diff2 = table.calc_diff(mr_calc, ref_calcs)

        assert diff1 == pytest.approx(0.00, abs=0.001)
        assert diff2 == pytest.approx(3.018, abs=0.001)


