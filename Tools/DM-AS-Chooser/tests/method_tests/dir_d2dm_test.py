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
Tests for DirD2DMTable.
"""
import pytest

from src.as_methods.dir_d2dm import DirD2DMTable
from src.utils import CONSTANTS
from tests.method_tests.mock_utils import FakeElectronicState, FakeMultiRefCalc, FakeSingleRefCalc, FakeVector


@pytest.fixture
def table():
    return DirD2DMTable(mr_calcs=[], reference_sources={})


class TestCallConvention:
    def test_calling_on_the_bare_class_raises_typeerror(self):
        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(0, 1.0, FakeVector(0, 0, 1.0)),
            FakeElectronicState(1, 3.0, FakeVector(0, 0, 3.0)),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.0 /CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 1.0 /CONSTANTS['DEBYE_PER_AU']))),
            FakeSingleRefCalc(FakeElectronicState(1, 3.0 /CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 3.0 /CONSTANTS['DEBYE_PER_AU']))),
        ]

        with pytest.raises(TypeError):
            DirD2DMTable.calc_diff(mr_calc, ref_calcs)

    def test_calling_on_an_instance_does_not_raise_typeerror(self, table):
        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(0, 1.0, FakeVector(0, 0, 1.0)),
            FakeElectronicState(1, 3.0, FakeVector(0, 0, 3.0)),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.0 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 1.0 /CONSTANTS['DEBYE_PER_AU']))),
            FakeSingleRefCalc(FakeElectronicState(1, 3.0 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 3.0 /CONSTANTS['DEBYE_PER_AU']))),
        ]

        try:
            table.calc_diff(mr_calc, ref_calcs)
        except TypeError:
            pytest.fail('calc_diff should accept (mr_calc, ref_calcs) when bound to an instance')
        except ValueError:
            pass


class TestCalcDiff:
    """
    Tests dir_gdm's calc_diff method. Should have identical logic to edm, except that it compares vectors
    rather than magnitudes.
    """

    def test_null_ground_cand(self, table):

        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(0, None, vector=None),
            FakeElectronicState(1, 5, FakeVector(0, 0, 6.0)),
            FakeElectronicState(2, 5, FakeVector(0, 0, 2.2)),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.0 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 1.0))),
            FakeSingleRefCalc(FakeElectronicState(1, 1.5 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 1.5))),
            FakeSingleRefCalc(FakeElectronicState(2, 3.0 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 3.0))),
        ]
        try:
            table.calc_diff(mr_calc, ref_calcs)
        except TypeError:
            pytest.fail('no ground-state entry found among candidate states')
        except ValueError:
            pass

    def test_null_excited_cand(self, table):

        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(0, 1.0, FakeVector(0, 0, 1.0)),
            FakeElectronicState(1, 6.0, None),
            FakeElectronicState(2, 2.2, None),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.0 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 1.0))),
            FakeSingleRefCalc(FakeElectronicState(1, 1.5 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 1.5))),
            FakeSingleRefCalc(FakeElectronicState(2, 3.0 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 3.0))),
        ]
        try:
            table.calc_diff(mr_calc, ref_calcs)
        except TypeError:
            pytest.fail('no excited-state dipole values found among candidate states')
        except ValueError:
            pass

    def test_null_ground_ref(self, table):

        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(0, 1.0, FakeVector(0, 0, 1.0)),
            FakeElectronicState(1, 6.0, FakeVector(0, 0, 6.0)),
            FakeElectronicState(2, 2.2, FakeVector(0, 0, 2.2)),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.0 / CONSTANTS['DEBYE_PER_AU'], None)),
            FakeSingleRefCalc(FakeElectronicState(1, 1.5 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 1.5))),
            FakeSingleRefCalc(FakeElectronicState(2, 3.0 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 3.0))),
        ]
        try:
            table.calc_diff(mr_calc, ref_calcs)
        except TypeError:
            pytest.fail('no ground-state entry found among reference states')
        except ValueError:
            pass

    def test_null_excited_ref(self, table):

        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(0, 1.0, FakeVector(0, 0, 1.0)),
            FakeElectronicState(1, 6.0, FakeVector(0, 0, 6.0)),
            FakeElectronicState(2, 2.2, FakeVector(0, 0, 2.2)),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.0 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 1.0))),
            FakeSingleRefCalc(FakeElectronicState(1, 1.5 / CONSTANTS['DEBYE_PER_AU'], None)),
            FakeSingleRefCalc(FakeElectronicState(2, 3.0 / CONSTANTS['DEBYE_PER_AU'], None)),
        ]
        try:
            table.calc_diff(mr_calc, ref_calcs)
        except TypeError:
            pytest.fail('no excited-state dipole values found among reference states')
        except ValueError:
            pass

    def test_matching_states_of_interest_gives_zero(self, table):
        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(0, 1.0, FakeVector(0, 0, 1.0)),
            FakeElectronicState(1, 1.5, FakeVector(0, 0, 1.5)),
            FakeElectronicState(2, 3.0, FakeVector(0, 0, 3.0)),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.0 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 1.0))),
            FakeSingleRefCalc(FakeElectronicState(1, 1.5 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 1.5))),
            FakeSingleRefCalc(FakeElectronicState(2, 3.0 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 3.0))),
        ]

        assert table.calc_diff(mr_calc, ref_calcs) == pytest.approx(0.0, abs=1e-6)

    def test_states_of_interest_need_not_share_state_number(self, table):
        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(0, 1.0, FakeVector(0, 0, 1.0)),
            FakeElectronicState(1, 6.0, FakeVector(0, 0, 6.0)),
            FakeElectronicState(2, 2.2, FakeVector(0, 0, 2.2)),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.0 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 1.0))),
            FakeSingleRefCalc(FakeElectronicState(1, 2.1 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 2.1))),
            FakeSingleRefCalc(FakeElectronicState(2, 6.0 / CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 6.0))),
        ]

        assert table.calc_diff(mr_calc, ref_calcs) == pytest.approx(0.0, abs=1e-4)
