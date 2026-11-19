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
Tests for orchestrator.py (file-discovery, reference-resolution, and symlink-creation helpers).
"""
import os

import pytest

from src.workflow_orchestrator import _matches_keyword, discover_files, create_excerr_symlinks, build_ref_calcs


class TestMatchesKeyword:
    """_matches_keyword: delimited-token matching, not unconstrained substring matching."""

    def test_matches_exact_token(self):
        assert _matches_keyword('td-b3lyp-01.orca', 'b3lyp')

    def test_does_not_match_substring_inside_longer_word(self):
        # "b3lyp" is a literal substring of "camb3lyp", but should NOT match
        assert not _matches_keyword('td-camb3lyp-01.orca', 'b3lyp')

    def test_m06_does_not_match_m062x(self):
        assert not _matches_keyword('td-m062x-01.orca', 'm06')

    def test_m06_does_not_match_m06l(self):
        assert not _matches_keyword('td-m06l-01.orca', 'm06')

    def test_pbe_does_not_match_pbe0(self):
        assert not _matches_keyword('td-pbe0-01.orca', 'pbe')

    def test_pbe_does_not_match_lcwhpbe(self):
        assert not _matches_keyword('td-lcwhpbe-01.orca', 'pbe')

    def test_case_insensitive(self):
        assert _matches_keyword('CASCI-06-06.LOG', 'casci')

    def test_casci_does_not_match_casscf(self):
        assert not _matches_keyword('casscf-06-06.log', 'casci')


class TestDiscoverFiles:
    """discover_files: groups files by keyword, filtered by extension."""

    @pytest.fixture
    def populated_dir(self, tmp_path):
        (tmp_path / 'casci-06-06.log').touch()
        (tmp_path / 'casci-06-07.log').touch()
        (tmp_path / 'casscf-08-07.csv').touch()
        (tmp_path / 'td-b3lyp-01.orca').touch()
        (tmp_path / 'td-b3lyp-02.out').touch()
        (tmp_path / 'td-camb3lyp-01.orca').touch()
        (tmp_path / 'ref-excs.list').touch()  # should never be picked up
        (tmp_path / 'notes.txt').touch()      # should never be picked up
        return tmp_path

    def test_groups_by_keyword(self, populated_dir):
        result = discover_files(str(populated_dir), ['casci', 'casscf'], valid_exts=('.log', '.csv'))
        assert set(result.keys()) == {'CASCI', 'CASSCF'}
        assert len(result['CASCI']) == 2
        assert len(result['CASSCF']) == 1

    def test_extension_filter_excludes_wrong_types(self, populated_dir):
        result = discover_files(str(populated_dir), ['casci', 'casscf'], valid_exts=('.log', '.csv'))
        all_matched = [f for files in result.values() for f in files]
        assert not any('ref-excs.list' in f for f in all_matched)
        assert not any('notes.txt' in f for f in all_matched)

    def test_reference_extension_filter(self, populated_dir):
        result = discover_files(str(populated_dir), ['b3lyp', 'camb3lyp'],
                                valid_exts=('.orca', '.out', '.csv'))
        assert set(result.keys()) == {'B3LYP', 'CAMB3LYP'}
        # b3lyp keyword does not pick up camb3lyp's files
        assert len(result['B3LYP']) == 2   # .orca and .out
        assert len(result['CAMB3LYP']) == 1

    def test_keyword_with_no_matches_is_omitted(self, populated_dir):
        result = discover_files(str(populated_dir), ['casci', 'nonexistent'], valid_exts=('.log',))
        assert 'NONEXISTENT' not in result

    def test_no_matches_returns_empty_dict(self, tmp_path):
        result = discover_files(str(tmp_path), ['casci'], valid_exts=('.log',))
        assert result == {}


class TestCreateExcerrSymlinks:
    """create_excerr_symlinks: winner -> .excerr.tsv, with graceful handling of missing data."""

    class FakeCalc:
        def __init__(self, num_electrons, num_orbitals, path):
            self.num_electrons = num_electrons
            self.num_orbitals = num_orbitals
            self.path = path

    @pytest.fixture
    def output_dir(self, tmp_path):
        extraction_dir = tmp_path / 'excited_state_extractions'
        extraction_dir.mkdir()
        (extraction_dir / 'casscf-06-06.excerr.tsv').write_text('fake\n')
        (extraction_dir / 'casscf-06-14.excerr.tsv').write_text('fake\n')
        # deliberately no .excerr.tsv for 06-08, to test graceful skip
        return tmp_path

    @pytest.fixture
    def mr_calcs(self):
        return [
            self.FakeCalc(6, 6, 'casscf-06-06.log'),
            self.FakeCalc(6, 14, 'casscf-06-14.log'),
            self.FakeCalc(6, 8, 'casscf-06-08.log'),
        ]

    def test_creates_symlink_for_winner_with_existing_excerr(self, output_dir, mr_calcs):
        winners = {'BLYP': '06-06'}
        create_excerr_symlinks('gdm-as', winners, mr_calcs, output_dir)

        link_path = output_dir / 'BLYP_gdm'
        assert link_path.is_symlink()
        assert link_path.resolve() == (output_dir / 'excited_state_extractions' / 'casscf-06-06.excerr.tsv').resolve()

    def test_symlink_target_is_absolute_to_prevent_broken_links(self, output_dir, mr_calcs):
        winners = {'BLYP': '06-06'}
        create_excerr_symlinks('gdm-as', winners, mr_calcs, output_dir)

        link_path = output_dir / 'BLYP_gdm'
        target_str = os.readlink(link_path)
        assert os.path.isabs(target_str)

    def test_skips_winner_with_no_excerr_file(self, output_dir, mr_calcs, caplog):
        winners = {'M06': '06-08'}  # 06-08 has no .excerr.tsv (see fixture)
        create_excerr_symlinks('gdm-as', winners, mr_calcs, output_dir)

        assert not (output_dir / 'M06_gdm').exists()

    def test_skips_none_winner(self, output_dir, mr_calcs):
        winners = {'PBE0': None}
        create_excerr_symlinks('gdm-as', winners, mr_calcs, output_dir)

        assert not (output_dir / 'PBE0_gdm').exists()

    def test_replaces_existing_symlink(self, output_dir, mr_calcs):
        winners = {'BLYP': '06-06'}
        create_excerr_symlinks('gdm-as', winners, mr_calcs, output_dir)
        # run again -- should not raise FileExistsError
        create_excerr_symlinks('gdm-as', winners, mr_calcs, output_dir)
        assert (output_dir / 'BLYP_gdm').is_symlink()

    def test_symlink_name_uses_method_without_as_suffix(self, output_dir, mr_calcs):
        winners = {'BLYP': '06-06'}
        create_excerr_symlinks('gdm-as', winners, mr_calcs, output_dir)
        assert (output_dir / 'BLYP_gdm').exists()
        assert not (output_dir / 'BLYP_gdm-as').exists()


class TestBuildRefCalcs:
    """
    build_ref_calcs: resolves raw reference-source paths into fully-pooled SingleRefCalc objects
    (unfiltered), mirroring build_mr_calcs' skip-and-warn-on-error pattern.
    """

    def test_pools_every_state_from_a_reference_file(self, tmp_path):
        ref_csv = tmp_path / 'pbe-01.csv'
        ref_csv.write_text('state,dipole\n0,1.98\n1,0.90\n')

        resolved = build_ref_calcs({'PBE': [str(ref_csv)]})

        assert len(resolved['PBE']) == 2
        assert {c.state.state for c in resolved['PBE']} == {0, 1}

    def test_pools_across_multiple_files_for_one_source(self, tmp_path):
        ref_csv1 = tmp_path / 'pbe-01.csv'
        ref_csv1.write_text('state,dipole\n0,1.98\n')
        ref_csv2 = tmp_path / 'pbe-02.csv'
        ref_csv2.write_text('state,dipole\n1,0.90\n')

        resolved = build_ref_calcs({'PBE': [str(ref_csv1), str(ref_csv2)]})

        assert len(resolved['PBE']) == 2
        assert {c.state.state for c in resolved['PBE']} == {0, 1}

    def test_skips_unreadable_file_with_a_warning_rather_than_raising(self, tmp_path):
        ref_csv = tmp_path / 'pbe-01.csv'
        ref_csv.write_text('state,dipole\n0,1.98\n')
        bad_path = tmp_path / 'nonexistent.csv'

        resolved = build_ref_calcs({'PBE': [str(ref_csv), str(bad_path)]})

        assert len(resolved['PBE']) == 1  # only the good file's state

    def test_multiple_sources_kept_independent(self, tmp_path):
        pbe_csv = tmp_path / 'pbe-01.csv'
        pbe_csv.write_text('state,dipole\n0,1.98\n')
        blyp_csv = tmp_path / 'blyp-01.csv'
        blyp_csv.write_text('state,dipole\n0,2.50\n')

        resolved = build_ref_calcs({'PBE': [str(pbe_csv)], 'BLYP': [str(blyp_csv)]})

        assert resolved['PBE'][0].state.dipole == pytest.approx(1.98)
        assert resolved['BLYP'][0].state.dipole == pytest.approx(2.50)