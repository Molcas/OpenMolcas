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
Tests for src/as_methods/base.py's shared functions: check_consistent_ground_states, filter_to_ground,
filter_to_excited.
"""
import pytest

from src.as_methods.base import check_consistent_ground_states, filter_to_ground, filter_to_excited


class FakeElectronicState:
    def __init__(self, state, dipole, vector=None):
        self.state = state
        self.dipole = dipole
        self.vector = vector


class FakeMultiRefCalc:
    def __init__(self, states, path='mr.log'):
        self.states = states
        self.path = path


class FakeSingleRefCalc:
    def __init__(self, state, path='ref.orca'):
        self.state = state
        self.path = path


class TestCheckConsistentGroundStates:
    def test_single_entry_never_raises(self):
        check_consistent_ground_states([FakeSingleRefCalc(FakeElectronicState(0, 1.98), 'file1')])

    def test_zero_entries_never_raises(self):
        check_consistent_ground_states([])

    def test_agreeing_entries_do_not_raise(self):
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 1.98), f'file{i}') for i in range(5)]
        check_consistent_ground_states(ref_calcs)

    def test_disagreeing_magnitude_raises(self):
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.98), 'file1'),
            FakeSingleRefCalc(FakeElectronicState(0, 9.0), 'file2'),
        ]
        with pytest.raises(ValueError, match='inconsistent ground-state dipole'):
            check_consistent_ground_states(ref_calcs)

    def test_same_magnitude_different_direction_raises(self):
        from types import SimpleNamespace
        vec1 = SimpleNamespace(x=0.0, y=0.0, z=-1.95)
        vec2 = SimpleNamespace(x=1.95, y=0.0, z=0.0)
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.95, vec1), 'file1'),
            FakeSingleRefCalc(FakeElectronicState(0, 1.95, vec2), 'file2'),
        ]
        with pytest.raises(ValueError, match='inconsistent ground-state vector'):
            check_consistent_ground_states(ref_calcs)

    def test_missing_vector_on_one_side_only_checks_magnitude(self):
        from types import SimpleNamespace
        vec1 = SimpleNamespace(x=0.0, y=0.0, z=-1.95)
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 1.95, vec1), 'file1'),
            FakeSingleRefCalc(FakeElectronicState(0, 1.95, None), 'file2'),
        ]
        check_consistent_ground_states(ref_calcs)  # should not raise -- vector check skipped safely


class TestFilterToGround:
    def test_narrows_candidate_to_ground_state_only(self):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 2.0), FakeElectronicState(1, 3.0)])
        mr_grounds, ref_grounds = filter_to_ground(mr_calc, [])
        assert [s.state for s in mr_grounds] == [0]

    def test_narrows_reference_to_ground_state_only(self):
        mr_calc = FakeMultiRefCalc([])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 2.0)),
                    FakeSingleRefCalc(FakeElectronicState(1, 3.0))]
        mr_grounds, ref_grounds = filter_to_ground(mr_calc, ref_calcs)
        assert len(ref_grounds) == 1
        assert ref_grounds[0].state.state == 0

    def test_collapses_multiple_consistent_ground_entries_to_one(self):
        mr_calc = FakeMultiRefCalc([])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 2.0), 'file1'),
            FakeSingleRefCalc(FakeElectronicState(0, 2.0), 'file2'),
            FakeSingleRefCalc(FakeElectronicState(1, 3.0), 'file1'),
        ]
        mr_grounds, ref_grounds = filter_to_ground(mr_calc, ref_calcs)
        assert len(ref_grounds) == 1

    def test_inconsistent_pooled_ground_entries_raise(self):
        mr_calc = FakeMultiRefCalc([])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 2.0), 'file1'),
            FakeSingleRefCalc(FakeElectronicState(0, 9.0), 'file2'),
        ]
        with pytest.raises(ValueError):
            filter_to_ground(mr_calc, ref_calcs)

    def test_does_not_mutate_original_ref_calcs(self):
        mr_calc = FakeMultiRefCalc([])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(0, 2.0), 'file1'),
            FakeSingleRefCalc(FakeElectronicState(0, 2.0), 'file2'),
        ]
        filter_to_ground(mr_calc, ref_calcs)
        assert len(ref_calcs) == 2, 'the original list must be untouched'


class TestFilterToExcited:
    def test_narrows_candidate_to_excited_states_only(self):
        mr_calc = FakeMultiRefCalc([FakeElectronicState(0, 2.0), FakeElectronicState(1, 3.0),
                                    FakeElectronicState(2, 4.0)])
        mr_excited, ref_excited = filter_to_excited(mr_calc, [])
        assert [s.state for s in mr_excited] == [1, 2]

    def test_narrows_reference_to_excited_states_only(self):
        mr_calc = FakeMultiRefCalc([])
        ref_calcs = [FakeSingleRefCalc(FakeElectronicState(0, 2.0)),
                    FakeSingleRefCalc(FakeElectronicState(1, 3.0))]
        mr_excited, ref_excited = filter_to_excited(mr_calc, ref_calcs)
        assert len(ref_excited) == 1
        assert ref_excited[0].state.state == 1

    def test_does_not_deduplicate_repeated_state_numbers(self):
        # documents a known, deliberately-unfixed gap: filter_to_excited narrows by state TYPE only,
        # not by state NUMBER -- two pooled files both claiming state=1 both pass through unchanged
        mr_calc = FakeMultiRefCalc([])
        ref_calcs = [
            FakeSingleRefCalc(FakeElectronicState(1, 0.5), 'file1'),
            FakeSingleRefCalc(FakeElectronicState(1, 5.2), 'file2'),
        ]
        mr_excited, ref_excited = filter_to_excited(mr_calc, ref_calcs)
        assert len(ref_excited) == 2