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
Tests for OldDirEDMTable.
"""
import pytest

from src.as_methods.old_dir_edm import OldDirEDMTable
from src.utils import CONSTANTS
from tests.method_tests.mock_utils import FakeElectronicState, FakeMultiRefCalc, FakeSingleRefCalc, FakeVector



@pytest.fixture
def table():
    return OldDirEDMTable(mr_calcs=[], reference_sources={})


class TestCallConvention:
    def test_calling_on_the_bare_class_raises_typeerror(self):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(1, 0.5, FakeVector(0, 0, 0.5))])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(1, 0.5 /CONSTANTS['DEBYE_PER_AU'],
                                                          FakeVector(0, 0, 0.5 /CONSTANTS['DEBYE_PER_AU'])))]

        with pytest.raises(TypeError):
            OldDirEDMTable.calc_diff(mr_calc, ref_calcs)

    def test_calling_on_an_instance_does_not_raise_typeerror(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(1, 0.5, FakeVector(0, 0, 0.5))])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(1, 0.5 /CONSTANTS['DEBYE_PER_AU'],
                                                          FakeVector(0, 0, 0.5 /CONSTANTS['DEBYE_PER_AU'])))]

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

    def test_null_excited_cand(self, table):
        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(1, 0.5, None),
            FakeElectronicState(2, 0.9, FakeVector(0, 0, 0.9)),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(1, 0.5 /CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 0.5))),
            FakeSingleRefCalc(FakeElectronicState(2, 0.9 /CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 0.9))),
        ]

        with pytest.raises(ValueError):
            table.calc_diff(mr_calc, ref_calcs)

    def test_null_excited_ref(self, table):
        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(1, 0.5, FakeVector(0, 0, 0.5)),
            FakeElectronicState(2, 0.9, FakeVector(0, 0, 0.9)),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(1, 0.5 /CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 0.5))),
            FakeSingleRefCalc(FakeElectronicState(2, 0.9 /CONSTANTS['DEBYE_PER_AU'], None)),
        ]

        with pytest.raises(ValueError):
            table.calc_diff(mr_calc, ref_calcs)

    def test_missing_same_state_passes(self, table):
        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(1, 0.5, None),
            FakeElectronicState(2, 0.9, FakeVector(0, 0, 0.9)),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(1, 0.5 /CONSTANTS['DEBYE_PER_AU'], None)),
            FakeSingleRefCalc(FakeElectronicState(2, 0.9 /CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 0.9))),
        ]

        assert table.calc_diff(mr_calc, ref_calcs) == pytest.approx(0.0)

    def test_matching_excited_state_vectors_gives_zero(self, table):
        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(1, 0.5, FakeVector(0, 0, 0.5)),
            FakeElectronicState(2, 0.9, FakeVector(0, 0, 0.9)),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(1, 0.5 /CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 0.5))),
            FakeSingleRefCalc(FakeElectronicState(2, 0.9 /CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 0.9))),
        ]

        assert table.calc_diff(mr_calc, ref_calcs) == pytest.approx(0.0, abs=1e-6)

    def test_returns_max_norm_difference_across_states(self, table):
        mr_calc = FakeMultiRefCalc([
            FakeElectronicState(1, 0.5, FakeVector(0, 0, 0.5)),
            FakeElectronicState(2, 0.9, FakeVector(0, 0, 0.9)),
        ])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(1, 0.6 /CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 0.6))),
            FakeSingleRefCalc(FakeElectronicState(2, 0.4 /CONSTANTS['DEBYE_PER_AU'], FakeVector(0, 0, 0.4))),
        ]

        assert table.calc_diff(mr_calc, ref_calcs) == pytest.approx(0.5, rel=1e-4)
