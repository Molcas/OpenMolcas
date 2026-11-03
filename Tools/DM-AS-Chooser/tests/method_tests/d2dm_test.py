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
Tests for D2DMTable.
"""
import pytest

from src.as_methods.d2dm import D2DMTable
from tests.method_tests.mock_utils import FakeElectronicState, FakeMultiRefCalc, FakeSingleRefCalc

@pytest.fixture
def table():
    return D2DMTable(mr_calcs=[], reference_sources={})


class TestCallConvention:
    def test_calling_on_the_bare_class_raises_typeerror(self):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 1.0), FakeElectronicState(1, 3.0)])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 1.0)),
                    FakeSingleRefCalc(FakeElectronicState(1, 3.0))]

        with pytest.raises(TypeError):
            D2DMTable.calc_diff(mr_calc, ref_calcs)

    def test_calling_on_an_instance_does_not_raise_typeerror(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 1.0), FakeElectronicState(1, 3.0)])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 1.0)),
                    FakeSingleRefCalc(FakeElectronicState(1, 3.0))]

        try:
            table.calc_diff(mr_calc, ref_calcs)
        except TypeError:
            pytest.fail('calc_diff should accept (mr_calc, ref_calcs) when bound to an instance')
        except ValueError:
            pass


class TestCalcDiff:
    """
    Test D2DM's calc_diff: both sides independently pick their excited state with the
    largest absolute dipole difference from their own ground state, then those two states-of-interest are
    compared directly.
    """

    def test_null_ground_cand(self, table):

        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(0, None), FakeElectronicState(1, 1.5), FakeElectronicState(2, 3.0),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.0)),
            FakeSingleRefCalc(FakeElectronicState(1, 1.5)),
            FakeSingleRefCalc(FakeElectronicState(2, 3.0)),
        ]
        try:
            table.calc_diff(mr_calc, ref_calcs)
        except TypeError:
            pytest.fail('no ground-state entry found among candidate states')
        except ValueError:
            pass

    def test_null_excited_cand(self, table):

        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(0, 1.0), FakeElectronicState(1, None), FakeElectronicState(2, None),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.0)),
            FakeSingleRefCalc(FakeElectronicState(1, 1.5)),
            FakeSingleRefCalc(FakeElectronicState(2, 3.0)),
        ]
        try:
            table.calc_diff(mr_calc, ref_calcs)
        except TypeError:
            pytest.fail('no excited-state dipole values found among candidate states')
        except ValueError:
            pass

    def test_null_ground_ref(self, table):

        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(0, 1.0), FakeElectronicState(1, 1.5), FakeElectronicState(2, 3.0),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, None)),
            FakeSingleRefCalc(FakeElectronicState(1, 1.5)),
            FakeSingleRefCalc(FakeElectronicState(2, 3.0)),
        ]
        try:
            table.calc_diff(mr_calc, ref_calcs)
        except TypeError:
            pytest.fail('no ground-state entry found among reference states')
        except ValueError:
            pass

    def test_null_excited_ref(self, table):

        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(0, 1.0), FakeElectronicState(1, 1.5), FakeElectronicState(2, 3.0),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.0)),
            FakeSingleRefCalc(FakeElectronicState(1, None)),
            FakeSingleRefCalc(FakeElectronicState(2, None)),
        ]
        try:
            table.calc_diff(mr_calc, ref_calcs)
        except TypeError:
            pytest.fail('no excited-state dipole values found among reference states')
        except ValueError:
            pass

    def test_no_ground_state_among_candidate_states_raises(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(1, 0.5), FakeElectronicState(2, 0.9)])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 1.98)),
                    FakeSingleRefCalc(FakeElectronicState(1, 0.5))]

        with pytest.raises(ValueError, match='candidate'):
            table.calc_diff(mr_calc, ref_calcs)

    def test_no_ground_state_among_ref_calcs_raises(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 1.98), FakeElectronicState(1, 0.5)])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(1, 0.5)),
                    FakeSingleRefCalc(FakeElectronicState(2, 0.9))]

        with pytest.raises(ValueError, match='reference'):
            table.calc_diff(mr_calc, ref_calcs)

    def test_no_excited_states_among_candidates_raises(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 1.98)])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 1.98)),
                    FakeSingleRefCalc(FakeElectronicState(1, 0.5))]

        with pytest.raises(ValueError, match='candidate'):
            table.calc_diff(mr_calc, ref_calcs)

    def test_selects_state_with_largest_dipole_difference(self, table):
        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(0, 1.0), FakeElectronicState(1, 1.5), FakeElectronicState(2, 3.0),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.0)),
            FakeSingleRefCalc(FakeElectronicState(1, 1.5)),
            FakeSingleRefCalc(FakeElectronicState(2, 3.0)),
        ]

        assert table.calc_diff(mr_calc, ref_calcs) == pytest.approx(0.0)

    def test_states_of_interest_dont_share_state_number(self, table):
        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(0, 1.0), FakeElectronicState(1, 5.0), FakeElectronicState(2, 1.2),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.0)),
            FakeSingleRefCalc(FakeElectronicState(1, 1.1)),
            FakeSingleRefCalc(FakeElectronicState(2, 5.0)),
        ]

        assert table.calc_diff(mr_calc, ref_calcs) == pytest.approx(0.0)