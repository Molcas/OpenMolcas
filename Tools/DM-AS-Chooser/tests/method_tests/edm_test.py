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
Tests for src/as_methods/edm.py's EDMTable.calc_diff.
"""
import pytest

from src.as_methods.edm import EDMTable

from tests.method_tests.mock_utils import FakeElectronicState, FakeMultiRefCalc, FakeSingleRefCalc

@pytest.fixture
def table():
    return EDMTable(mr_calcs=[], reference_sources={})


class TestCallConvention:
    def test_calling_on_the_bare_class_raises_typeerror(self):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(1, 0.5)])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(1, 0.5))]

        with pytest.raises(TypeError):
            EDMTable.calc_diff(mr_calc, ref_calcs)

    def test_calling_on_an_instance_does_not_raise_typeerror(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(1, 0.5)])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(1, 0.5))]

        try:
            table.calc_diff(mr_calc, ref_calcs)
        except TypeError:
            pytest.fail('calc_diff should accept (mr_calc, ref_calcs) when bound to an instance')
        except ValueError:
            pass


class TestCalcDiff:
    def test_matching_excited_state_dipoles_gives_zero(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(1, 0.5), FakeElectronicState(2, 0.9)])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(1, 0.5)),
                    FakeSingleRefCalc(FakeElectronicState(2, 0.9))]

        assert table.calc_diff(mr_calc, ref_calcs) == pytest.approx(0.0)

    def test_returns_max_absolute_difference_across_states(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(1, 0.5), FakeElectronicState(2, 0.9)])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(1, 0.6)),
                    FakeSingleRefCalc(FakeElectronicState(2, 0.4))]

        assert table.calc_diff(mr_calc, ref_calcs) == pytest.approx(0.5)

    def test_state_set_mismatch_raises(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(1, 0.5), FakeElectronicState(2, 0.9)])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(1, 0.5)),
                    FakeSingleRefCalc(FakeElectronicState(3, 0.9))]

        with pytest.raises(ValueError, match='state mismatch'):
            table.calc_diff(mr_calc, ref_calcs)

    def test_calc_diff_has_ground_state_awareness(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 2.0), FakeElectronicState(1, 0.5)])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 9.0)),
                    FakeSingleRefCalc(FakeElectronicState(1, 0.5))]

        result = table.calc_diff(mr_calc, ref_calcs)
        assert result == pytest.approx(0.0)