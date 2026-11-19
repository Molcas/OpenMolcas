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
Tests for src/as_methods/old_dir_gdm.py's OldDirGDMTable.calc_diff.
"""
import pytest

from src.as_methods.old_dir_gdm import OldDirGDMTable

from tests.method_tests.mock_utils import FakeElectronicState, FakeMultiRefCalc, FakeSingleRefCalc, FakeVector

@pytest.fixture
def table():
    return OldDirGDMTable(mr_calcs=[], reference_sources={})


class TestCallConvention:
    def test_calling_on_the_bare_class_raises_typeerror(self):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 1.95, FakeVector(0, 0, -1.95))])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 1.95, FakeVector(0, 0, -1.95)))]

        with pytest.raises(TypeError):
            OldDirGDMTable.calc_diff(mr_calc, ref_calcs)

    def test_calling_on_an_instance_does_not_raise_typeerror(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 1.95, FakeVector(0, 0, -1.95))])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 1.95, FakeVector(0, 0, -1.95)))]

        try:
            table.calc_diff(mr_calc, ref_calcs)
        except TypeError:
            pytest.fail('calc_diff should accept (mr_calc, ref_calcs) when bound to an instance')
        except ValueError:
            pass


class TestCalcDiff:
    def test_matching_ground_state_vectors_gives_zero(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 1.95, FakeVector(0, 0, -1.95))])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 1.95, FakeVector(0, 0, -1.95)))]

        assert table.calc_diff(mr_calc, ref_calcs) == pytest.approx(0.0, abs=1e-6)

    def test_differing_vectors_gives_norm_of_difference(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 1.95, FakeVector(0, 0, -1.95))])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 0.767, FakeVector(0, 0, -0.767)))]

        expected = abs(-1.95 - (-0.767))
        assert table.calc_diff(mr_calc, ref_calcs) == pytest.approx(expected, rel=1e-4)


class TestFilterToGroundSharedWithGdm:
    def test_collapses_multiple_consistent_ground_entries_pooled_across_files(self, table):
        # confirms OldDirGDMTable actually gets this behavior via the shared
        # filter_to_ground function, not a duplicated copy of its own
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 1.95, FakeVector(0, 0, -1.95))])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.95, FakeVector(0, 0, -1.95)), 'file1'),
            FakeSingleRefCalc(FakeElectronicState(0, 1.95, FakeVector(0, 0, -1.95)), 'file2'),
        ]

        assert table.calc_diff(mr_calc, ref_calcs) == pytest.approx(0.0, abs=1e-6)

    def test_inconsistent_pooled_ground_entries_raise(self, table):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 1.95, FakeVector(0, 0, -1.95))])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.95, FakeVector(0, 0, -1.95)), 'file1'),
            FakeSingleRefCalc(FakeElectronicState(0, 5.0, FakeVector(0, 0, -5.0)), 'file2'),
        ]

        with pytest.raises(ValueError, match='inconsistent'):
            table.calc_diff(mr_calc, ref_calcs)
