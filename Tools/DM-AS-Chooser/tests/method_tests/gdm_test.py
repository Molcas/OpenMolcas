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
Tests for src/as_methods/gdm.py's GDMTable.calc_diff.
"""
import pytest

from src.as_methods.gdm import GDMTable

from tests.method_tests.mock_utils import FakeElectronicState, FakeMultiRefCalc, FakeSingleRefCalc

@pytest.fixture
def table():
    return GDMTable(mr_calcs=[], reference_sources={})


class TestCallConvention:
    def test_calling_on_the_bare_class_raises_typeerror(self):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 2.0)])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 2.0))]

        with pytest.raises(TypeError):
            GDMTable.calc_diff(mr_calc, ref_calcs)

    def test_calling_on_an_instance_does_not_raise_typeerror(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 2.0)])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 2.0))]

        try:
            table.calc_diff(mr_calc, ref_calcs)
        except TypeError:
            pytest.fail('calc_diff should accept (mr_calc, ref_calcs) when bound to an instance')
        except ValueError:
            pass


class TestCalcDiff:
    def test_matching_ground_state_dipoles_gives_zero(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 2.0), FakeElectronicState(1, 3.0)])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 2.0))]

        assert table.calc_diff(mr_calc, ref_calcs) == pytest.approx(0.0)

    def test_differing_ground_state_dipoles_gives_absolute_difference(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 2.0), FakeElectronicState(1, 3.0)])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 2.5))]

        assert table.calc_diff(mr_calc, ref_calcs) == pytest.approx(0.5)

    def test_zero_ref_calcs_raises(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 2.0)])

        with pytest.raises(ValueError, match='exactly one'):
            table.calc_diff(mr_calc, [])

    def test_inconsistent_pooled_ground_entries_raise(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 2.0)])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 2.0), 'file1'),
            FakeSingleRefCalc(FakeElectronicState(0, 9.0), 'file2'),
        ]

        with pytest.raises(ValueError, match='inconsistent ground-state dipole'):
            table.calc_diff(mr_calc, ref_calcs)

    def test_no_ground_state_among_candidate_states_raises(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(1, 3.0)])  # no state==0
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 2.0))]

        with pytest.raises(ValueError, match='expected exactly one ground-state entry, got 0 in mr.log'):
            table.calc_diff(mr_calc, ref_calcs)