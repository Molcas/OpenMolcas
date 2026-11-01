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
Tests for src/generate_molcas_files.py's build_timing_input, unique_winners, and generate_timing_inputs.
"""
import pytest

from src.generate_molcas_files import build_comparison_input, unique_winners, generate_comparison_inputs


class FakeArgs:
    coord = "starting.xyz"
    basis = "ANO-RCC-VTZP"
    integrals = "$CurrDir/2e-integrals/2e-integrals.OrdInt"
    fileorb = "$CurrDir/uhf/uhf.UnaOrb"
    nroots = 6
    charge = 0
    spin = 1
    cimx = 200
    caspt2_maxiter = 10000
    ksdft = "T:PBE"
    fix_caspt2_grid_loop_var = False


@pytest.fixture
def args():
    return FakeArgs()


class TestBuildTimingInput:
    def test_substitutes_nael_and_norb(self, args):
        result = build_comparison_input(6, 8, args)
        assert 'nActEl = 6' in result
        assert 'Ras2 = 8' in result

    def test_different_active_space_gives_different_content(self, args):
        result1 = build_comparison_input(6, 6, args)
        result2 = build_comparison_input(10, 12, args)
        assert result1 != result2

    def test_does_not_contain_cionly(self, args):
        # this is the key distinguishing feature from the selection-phase
        # CASCI input, which uses CIOnly -- this is a full CASSCF calculation
        result = build_comparison_input(6, 8, args)
        assert 'CIOnly' not in result

    def test_contains_all_expected_sections(self, args):
        result = build_comparison_input(6, 8, args)
        for section in ('&GATEWAY', '&SEWARD', '&RASSCF', '&GRID_IT', '&CASPT2', '&MCPDFT', '&RASSI'):
            assert section in result


    def test_fix_flag_uses_correct_loop_var(self, args):
        args.fix_caspt2_grid_loop_var = True
        result = build_comparison_input(6, 8, args)
        assert 'FileOrb = $Project.Pt2Orb.$B' in result
        assert 'FileOrb = $Project.Pt2Orb.$A' not in result

    def test_ciroot_uses_nroots_for_both_positions(self, args):
        args.nroots = 4
        result = build_comparison_input(6, 8, args)
        assert 'CIRoot = 4 4 1' in result

    def test_substitutes_ksdft_functional(self, args):
        args.ksdft = 'T:BLYP'
        result = build_comparison_input(6, 8, args)
        assert 'KSDFT = T:BLYP' in result


class TestUniqueWinners:
    def test_deduplicates_across_methods_and_sources(self):
        method_results = {
            'gdm-as': {'BLYP': '06-06', 'PBE': '06-06'},
            'edm-as': {'BLYP': '06-07', 'PBE': '06-06'},
        }
        assert unique_winners(method_results) == ['06-06', '06-07']

    def test_skips_none_winners(self):
        method_results = {'gdm-as': {'BLYP': '06-06', 'PBE': None}}
        assert unique_winners(method_results) == ['06-06']

    def test_empty_results_gives_empty_list(self):
        assert unique_winners({}) == []

    def test_returns_sorted(self):
        method_results = {'gdm-as': {'BLYP': '10-10', 'PBE': '06-06', 'M06': '08-08'}}
        assert unique_winners(method_results) == ['06-06', '08-08', '10-10']
