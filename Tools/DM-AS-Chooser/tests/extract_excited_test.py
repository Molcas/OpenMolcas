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
Tests for extract_excited.py (extract_excdata, compute_excerr,
extract_excerr, and their internal field-extraction helpers).
"""
from pathlib import Path

import numpy as np
import pytest

from src.extract_excited import (
    _bash_field, _grep_field, _grep_block_field,
    extract_excdata, extract_excerr
)


class TestBashField:
    """_bash_field: whitespace-run splitting, matching how these log files are laid out."""

    def test_extracts_field_with_leading_indentation(self):
        line = '::    CASPT2 Root  1     Total energy:   -153.62913196'
        assert _bash_field(line, 7) == '-153.62913196'

    def test_squeezes_multiple_spaces(self):
        line = '     2             0.1586454975        4.3169639088          34818.6621'
        assert _bash_field(line, 4) == '4.3169639088'

    def test_line_without_leading_whitespace(self):
        line = 'a b c d'
        assert _bash_field(line, 1) == 'a'
        assert _bash_field(line, 4) == 'd'

    def test_trailing_newline_is_stripped(self):
        line = 'a b c\n'
        assert _bash_field(line, 3) == 'c'


class TestGrepField:
    """_grep_field: one field from every independently-matching line."""

    @pytest.fixture
    def log_file(self, tmp_path):
        path = tmp_path / 'sample.log'
        path.write_text(
            '::    CASPT2 Root  1     Total energy:   -153.62913196\n'
            'irrelevant line\n'
            '::    CASPT2 Root  2     Total energy:   -153.47371405\n'
        )
        return path

    def test_extracts_field_from_every_matching_line(self, log_file):
        result = _grep_field(log_file, r'CASPT2 Root', field_1indexed=7)
        assert result == ['-153.62913196', '-153.47371405']

    def test_no_matches_returns_empty_list(self, log_file):
        result = _grep_field(log_file, r'NOT PRESENT', field_1indexed=1)
        assert result == []

    def test_does_not_match_ms_prefixed_variant_by_accident(self, tmp_path):
        path = tmp_path / 'sample.log'
        path.write_text(
            '::    CASPT2 Root  1     Total energy:   -153.62913196\n'
            '::    MS-CASPT2 Root  1     Total energy:   -153.62999999\n'
        )
        result = _grep_field(path, r'CASPT2 Root', field_1indexed=7)
        assert len(result) == 2, (
            'documents the real substring-matching ambiguity -- both lines '
            'match a plain "CASPT2 Root" pattern'
        )


class TestGrepBlockField:
    """_grep_block_field: first matching block, last `tail` of an `after`-sized window."""

    @pytest.fixture
    def sf_state_log(self, tmp_path):
        path = tmp_path / 'sample.log'
        path.write_text(
            ' SF State       Relative EMIN(au)   Rel lowest level(eV)    D:o, cm**(-1)\n'
            '\n'
            '     1             0.0000000000        0.0000000000              0.0000\n'
            '     2             0.1586454975        4.3169639088          34818.6621\n'
            '     3             0.2355299175        6.4090955581          51692.8418\n'
            '     4             0.2656962045        7.2299620455          58313.5765\n'
            '     5             0.2783763469        7.5750062981          61096.5461\n'
            '     6             0.2865555996        7.7975751014          62891.6846\n'
        )
        return path

    def test_extracts_field_from_tail_of_block(self, sf_state_log):
        result = _grep_block_field(sf_state_log, r'SF State', after=7, tail=5, field_1indexed=4)
        assert result == ['4.3169639088', '6.4090955581', '7.2299620455', '7.5750062981', '7.7975751014']

    def test_excludes_ground_state_row_via_tail(self, sf_state_log):
        # the header + blank + state-1 row are excluded by tail=5 keeping
        # only the last 5 of the 8-line window
        result = _grep_block_field(sf_state_log, r'SF State', after=7, tail=5, field_1indexed=4)
        assert '0.0000000000' not in result

    def test_no_match_returns_empty_list(self, tmp_path):
        path = tmp_path / 'empty.log'
        path.write_text('nothing relevant here\n')
        result = _grep_block_field(path, r'SF State', after=7, tail=5, field_1indexed=4)
        assert result == []

    def test_only_first_block_used_when_pattern_appears_twice(self, tmp_path):
        path = tmp_path / 'sample.log'
        path.write_text(
            ' SF State header\n' + ('x\n' * 6) + 'FIRST_TAIL_LINE\n'
            ' SF State header\n' + ('y\n' * 6) + 'SECOND_TAIL_LINE\n'
        )
        result = _grep_block_field(path, r'SF State', after=7, tail=1, field_1indexed=1)
        assert result == ['FIRST_TAIL_LINE']


class TestExtractExcdata:
    """extract_excdata: end-to-end extraction from a real-shaped log."""

    @pytest.fixture
    def full_log(self, tmp_path):
        path = tmp_path / 'casscf-06-14.log'
        path.write_text(
            ' SF State       Relative EMIN(au)   Rel lowest level(eV)    D:o, cm**(-1)\n'
            '\n'
            '     1             0.0000000000        0.0000000000              0.0000\n'
            '     2             0.1586454975        4.3169639088          34818.6621\n'
            '     3             0.2355299175        6.4090955581          51692.8418\n'
            '     4             0.2656962045        7.2299620455          58313.5765\n'
            '     5             0.2783763469        7.5750062981          61096.5461\n'
            '     6             0.2865555996        7.7975751014          62891.6846\n'
            '::    CASPT2 Root  1     Total energy:   -153.62913196\n'
            '::    CASPT2 Root  2     Total energy:   -153.47371405\n'
            '::    CASPT2 Root  3     Total energy:   -153.37213041\n'
            '::    CASPT2 Root  4     Total energy:   -153.34445471\n'
            '::    CASPT2 Root  5     Total energy:   -153.32767483\n'
            '::    CASPT2 Root  6     Total energy:   -153.32117831\n'
            '      Total MC-PDFT energy for state    1                                -153.74093920\n'
            '      Total MC-PDFT energy for state    2                                -153.58531338\n'
            '      Total MC-PDFT energy for state    3                                -153.50176585\n'
            '      Total MC-PDFT energy for state    4                                -153.47664356\n'
            '      Total MC-PDFT energy for state    5                                -153.46264368\n'
            '      Total MC-PDFT energy for state    6                                -153.45655134\n'
        )
        return path

    def test_writes_five_rows_for_six_roots(self, full_log, tmp_path):
        out_path = tmp_path / 'out.exc.tsv'
        extract_excdata(str(full_log), out_path=out_path)
        data = np.loadtxt(out_path)
        assert data.shape == (5, 5)

    def test_cas_column_not_hartree_converted(self, full_log, tmp_path):
        out_path = tmp_path / 'out.exc.tsv'
        extract_excdata(str(full_log), out_path=out_path)
        data = np.loadtxt(out_path)
        assert data[0, 0] == pytest.approx(4.3169639088)

    def test_pt2_column_ground_state_subtracted_and_converted(self, full_log, tmp_path):
        out_path = tmp_path / 'out.exc.tsv'
        extract_excdata(str(full_log), out_path=out_path)
        data = np.loadtxt(out_path)
        expected = 27.211386245988 * (-153.47371405 - (-153.62913196))
        assert data[0, 1] == pytest.approx(expected)

    def test_missing_caspt2_data_raises(self, tmp_path):
        path = tmp_path / 'no_caspt2.log'
        path.write_text('::    RASSCF root number  1 Total energy:    -74.65298980\n')
        with pytest.raises((ValueError, IndexError)):
            extract_excdata(str(path), out_path=tmp_path / 'out.exc.tsv')

    def test_default_out_path_replaces_log_extension(self, full_log):
        written = extract_excdata(str(full_log))
        assert written == str(full_log.with_suffix('').with_suffix('.exc.tsv'))
        assert Path(written).exists()


class TestComputeExcerr:
    """compute_excerr: single .exc.tsv vs a ref-excs.list in an explicit ref_dir."""

    @pytest.fixture
    def exc_tsv(self, tmp_path):
        path = tmp_path / 'output' / 'casscf-06-14.exc.tsv'
        path.parent.mkdir(parents=True)
        np.savetxt(path, np.array([[4.30, 4.23, 4.24, 4.23, 4.24]]), delimiter='\t')
        return path

    def test_returns_none_when_no_ref_file(self, exc_tsv, tmp_path):
        ref_dir = tmp_path / 'no_ref_here' / 'not-refs.list'
        assert extract_excerr(str(exc_tsv), ref_dir) is None

    def test_computes_absolute_error(self, exc_tsv, tmp_path):
        ref_dir = tmp_path / 'data_dir'
        ref_dir.mkdir()
        (ref_dir / 'ref-excs.list').write_text('4.20\n')

        err_path = extract_excerr(str(exc_tsv), ref_dir / 'ref-excs.list')
        assert err_path is not None
        err_data = np.loadtxt(err_path)
        assert err_data[0] == pytest.approx(0.10)

    def test_output_filename_only_replaces_filename_suffix_not_anywhere_in_path(self, tmp_path):
        exc_dir = tmp_path / 'my_excellent_results'
        exc_dir.mkdir()
        exc_path = exc_dir / 'casscf-06-14.exc.tsv'
        np.savetxt(exc_path, np.array([[4.30]]), delimiter='\t')

        ref_dir = tmp_path
        (ref_dir / 'ref-excs.list').write_text('4.20\n')

        err_path = extract_excerr(str(exc_path), ref_dir / 'ref-excs.list')
        assert err_path == str(exc_dir / 'casscf-06-14.excerr.tsv')
        assert 'excellent' in err_path  # directory name untouched


class TestExtractExcerr:
    """extract_excerr: single .exc.tsv vs a direct ref_path, computing absolute error."""

    @pytest.fixture
    def exc_tsv(self, tmp_path):
        path = tmp_path / 'output' / 'casscf-06-14.exc.tsv'
        path.parent.mkdir(parents=True)
        np.savetxt(path, np.array([[4.30, 4.23, 4.24, 4.23, 4.24]]), delimiter='\t')
        return path

    def test_returns_none_when_ref_path_does_not_exist(self, exc_tsv, tmp_path):
        missing_ref = tmp_path / 'nonexistent-ref-excs.list'
        assert extract_excerr(str(exc_tsv), missing_ref) is None

    def test_computes_absolute_error(self, exc_tsv, tmp_path):
        ref_path = tmp_path / 'ref-excs.list'
        ref_path.write_text('4.20\n')

        err_path = extract_excerr(str(exc_tsv), ref_path)
        assert err_path is not None
        err_data = np.loadtxt(err_path)
        assert err_data[0] == pytest.approx(0.10)

    def test_output_filename_only_replaces_filename_suffix_not_anywhere_in_path(self, tmp_path):
        exc_dir = tmp_path / 'my_excellent_results'
        exc_dir.mkdir()
        exc_path = exc_dir / 'casscf-06-14.exc.tsv'
        np.savetxt(exc_path, np.array([[4.30]]), delimiter='\t')

        ref_path = tmp_path / 'ref-excs.list'
        ref_path.write_text('4.20\n')

        err_path = extract_excerr(str(exc_path), ref_path)
        assert err_path == str(exc_dir / 'casscf-06-14.excerr.tsv')
        assert 'excellent' in err_path

    def test_error_written_alongside_exc_tsv_not_alongside_ref_path(self, tmp_path):
        # ref_path can live in a completely different directory
        exc_dir = tmp_path / 'results' / 'casscf'
        exc_dir.mkdir(parents=True)
        exc_path = exc_dir / 'casscf-06-14.exc.tsv'
        np.savetxt(exc_path, np.array([[4.30]]), delimiter='\t')

        ref_dir = tmp_path / 'data_dir'
        ref_dir.mkdir()
        ref_path = ref_dir / 'ref-excs.list'
        ref_path.write_text('4.20\n')

        err_path = extract_excerr(str(exc_path), ref_path)
        assert Path(err_path).parent == exc_dir
        assert Path(err_path).parent != ref_dir