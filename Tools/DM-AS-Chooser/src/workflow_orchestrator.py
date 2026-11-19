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


import csv
import os
import re
from pathlib import Path
from typing import Dict, List

from src.generate_molcas_files import generate_comparison_inputs
from src.as_methods.dir_d2dm import DirD2DMTable
from src.as_methods.dir_edm import DirEDMTable
from src.as_methods.old_dir_edm import OldDirEDMTable
from src.as_methods.dir_gdm import DirGDMTable
from src.as_methods.old_dir_gdm import OldDirGDMTable
from src.as_methods.edm import EDMTable
from src.as_methods.gdm import GDMTable
from src.as_methods.d2dm import D2DMTable
from src.calc_base import MultiRefCalc, SingleRefCalc
from src.choose import pick_winners
from src.extract_excited import extract_excdata, extract_excerr
from src.utils import (GDM_AS, EDM_AS, DipoleNotFoundError, logger, D2DM_AS, DIR_GDM_AS, DIR_EDM_AS, DIR_D2DM_AS, OLD_DIR_GDM_AS, OLD_DIR_EDM_AS,
                       LOG_EXT, CSV_EXT, VALID_ORCA_EXTS, EXTRACTABLE_KEYWORDS)

# Output filename of final active space selections
SELECTIONS_FILE_NAME = 'active_space_selections.tsv'

# name -> Table class.
PROTOCOL_SPECS = {
    GDM_AS: GDMTable,
    EDM_AS: EDMTable,
    D2DM_AS: D2DMTable,
    DIR_GDM_AS: DirGDMTable,
    DIR_EDM_AS: DirEDMTable,
    DIR_D2DM_AS: DirD2DMTable,
    OLD_DIR_GDM_AS: OldDirGDMTable,
    OLD_DIR_EDM_AS: OldDirEDMTable,
}


def _matches_keyword(filename: str, keyword: str) -> bool:
    """
    Checks whether `keyword` appears in `filename` as a delimited token (bounded by the start/end of the
    name or a non-alphanumeric character/delimiter on each side).

    :param filename: the filename to check.
    :param keyword: the keyword to search for.

    :returns: True if `keyword` appears in `filename` as a delimited token, False otherwise.
    """
    pattern = r'(?<![a-zA-Z0-9])' + re.escape(keyword) + r'(?![a-zA-Z0-9])'
    return re.search(pattern, filename, re.IGNORECASE) is not None


def discover_files(data_dir: str, keywords: List[str], valid_exts: tuple) -> Dict[str, List[str]]:
    """
    For each keyword, finds every file in `data_dir` (limited to `valid_exts`) whose name contains that
    keyword as a delimited token, grouped by keyword.

    :param data_dir: the directory to search.
    :param keywords: keyword(s) identifying files by name.
    :param valid_exts: which file extensions to consider; anything else is ignored.

    :returns: mapping of keyword (uppercased) -> list of matching file paths. A keyword with no matches
        is omitted from the result.
    """
    data_dir = Path(data_dir)
    all_files = sorted(
        str(p) for p in data_dir.iterdir()
        if p.is_file() and p.suffix in valid_exts
    )

    grouped = {}
    for kw in keywords:
        matches = [f for f in all_files if _matches_keyword(Path(f).name, kw)]
        if matches:
            grouped[kw.upper()] = matches

    return grouped


def build_mr_calcs(paths: List[str]) -> List[MultiRefCalc]:
    """
    Constructs a MultiRefCalc for each path, warning and skipping any file that fails to parse or validate.

    :param paths: a list of paths to multiconfiguration log file to parse into MultiRefCalcs.

    :returns: a list of MultiRefCalcs objects, one for each path.
    """
    mr_calcs = []
    for path in paths:
        try:
            mr_calcs.append(MultiRefCalc.infer(path))
        except (ValueError, DipoleNotFoundError) as e:
            logger.warning(f'{path}: {e}; removing from analysis')

    return mr_calcs

def build_ref_calcs(ref_sources: Dict[str, List[str]]) -> Dict[str, List[SingleRefCalc]]:
    """
    Resolves reference-source file paths into SingleRefCalc objects.

    :param ref_sources: mapping of reference method name -> list of raw file paths.
                        ex: {"PBE": ["pbe-file1", pbe_file2" ...], "BLYP": ...}

    :returns: mapping of reference method name -> list of resolved SingleRefCalc objects.
              ex: {"PBE": [SingleRefCalc(state=1, ...), ...], "BLYP": ...}
    """
    resolved = {}
    for name, paths in ref_sources.items():
        calcs = []
        for p in paths:
            try:
                calcs.extend(SingleRefCalc.all_from_file(p))
            except (ValueError, DipoleNotFoundError) as e:
                logger.warning(f'{p}: {e}; skipping')
        resolved[name] = calcs
    return resolved


def update_selections_file(method: str, winners: dict, output_dir: Path):
    """
    Records this method's winning active space per reference source
    in a shared selections file, e.g.:
        gdm    BLYP    06-06
        edm    BLYP    06-06

    :param method: the excited-state single reference method name be recorded.
    :param winners: the winning active space per reference source.
    :param output_dir: the path to a shared selections file.
    """
    label = method[:-3] if method.endswith('-as') else method

    path = output_dir / SELECTIONS_FILE_NAME

    rows = []
    if os.path.exists(path):
        with open(path, 'r', newline='') as f:
            reader = csv.reader(f, delimiter='\t')
            rows = [row for row in reader if row and row[0] != label]

    for source, winner in winners.items():
        rows.append([label, source, winner if winner is not None else ''])

    with open(path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(rows)


def create_excerr_symlinks(method: str, winners: dict, mr_calcs: List[MultiRefCalc], output_dir: Path):
    """
    For each reference source's winning active space, creates a symlink named after the reference source
    and protocol, pointing to that active space's .excerr.tsv file.

    PURPOSE: target of the symlink is the winning active space.

    :param method: the TD_DFT method name, used in the symlink name.
    :param winners: winning active space label per reference source.
    :param mr_calcs: the candidate MultiRefCalc objects, used to map a winning label back to its source
        log file's basename.
    :param output_dir: the directory containing both the .excerr.tsv files and where the symlinks are
        created.
    """
    label = method[:-3] if method.endswith('-as') else method

    label_to_stem = {
        f'{c.num_electrons:02d}-{c.num_orbitals:02d}': Path(c.path).stem
        for c in mr_calcs
    }

    extraction_dir = output_dir / 'excited_state_extractions'

    for source, winner in winners.items():
        if winner is None:
            continue

        stem = label_to_stem.get(winner)
        if stem is None:
            logger.warning(f'{source}/{label}: winning label {winner} not found among candidates; skipping symlink')
            continue

        target = extraction_dir / f'{stem}.excerr.tsv'
        if not target.exists():
            logger.warning(f'{source}/{label}: {target} does not exist; skipping symlink')
            continue

        link_path = output_dir / f'{source}_{label}'
        if link_path.exists() or link_path.is_symlink():
            link_path.unlink()
        link_path.symlink_to(target)


def run_method(method: str, parser, mr_calcs: List[MultiRefCalc],
               ref_calcs: Dict[str, List[SingleRefCalc]], output_dir: Path) -> dict:
    """
    Runs one protocol end-to-end against the shared multiconfigurational candidate active spaces and reference-source
    files, returning its winners dict.

    Writes a .tsv file for each protocol and records winning active spaces per protocol and method.

    :param method: the TD_DFT method name being run.
    :param parser: the argument parser object.
    :param mr_calcs: a list of MultiRefCalc objects, one for each candidate active space.
    :param ref_calcs: the map of reference methods to their objects {"PBE": [SingleRefCalc(state=1, ...), ...], "BLYP": ...}
    :param output_dir: the output directory to write to.

    :returns: a dictionary of winning active spaces per protocol and method.
    """
    table_cls = PROTOCOL_SPECS[method]

    try:
        table = table_cls(mr_calcs, ref_calcs)
    except ValueError as e:
        parser.error(str(e))
        return {}

    columns, rows = table.build()
    winners = pick_winners(columns, rows)

    update_selections_file(method, winners, output_dir)

    tsv_path = table.write_tsv(output_dir)
    logger.info(f'wrote table to {tsv_path}')

    return winners


def run_workflow(opts, parser) -> dict:
    """
    Runs the full active-space-selection workflow for the given parsed CLI options:
        1. Discovers active space candidate and reference files
        2. Builds and sorts MultiRefCalc objects per -m keyword group
        3. Optionally runs excitation-energy extraction/benchmarking
        4. Runs every requested protocol per keyword group and generates a .tsv file for each keyword (protocol)
        5. Records winners for each protocol and (if extraction ran) symlinks.

    :param opts: parsed CLI options (see main.get_parser).
    :param parser: the argument parser, used to report validation failures.

    :returns: mapping of keyword -> {method: winners}.
    """
    if opts.output_dir:
        base_output_dir = Path(opts.output_dir)
    else:
        base_output_dir = Path('results').resolve()

    base_output_dir.mkdir(parents=True, exist_ok=True)

    # Source the active space candidates
    mr_groups = discover_files(opts.data_dir, opts.mr_keywords, valid_exts=(LOG_EXT, CSV_EXT))
    if not mr_groups:
        parser.error(f'no *.log files matching {opts.mr_keywords} found in {opts.data_dir}')

    # Source the reference files
    ref_sources = discover_files(opts.data_dir, opts.ref_keywords, valid_exts=VALID_ORCA_EXTS + (CSV_EXT,))
    if not ref_sources:
        parser.error(f'no files matching {opts.ref_keywords} found in {opts.data_dir}')

    all_results = {}
    for keyword, mr_files in mr_groups.items():

        # Build the MultiRefCalc objects
        mr_calcs = build_mr_calcs(mr_files)
        ref_calcs = build_ref_calcs(ref_sources)
        if not mr_calcs:
            logger.warning(f'{keyword}: no valid multi reference calcs found; skipping')
            continue

        # Ensure .tsv file ordering
        mr_calcs.sort(key=lambda c: (c.num_electrons, c.num_orbitals))

        output_dir = base_output_dir / keyword.lower()
        output_dir.mkdir(parents=True, exist_ok=True)

        # Run excited-state data extraction (if requested and allowed)
        do_extract = opts.extract_excdata and keyword in EXTRACTABLE_KEYWORDS
        if do_extract:
            extraction_dir = output_dir / 'excited_state_extractions'
            extraction_dir.mkdir(parents=True, exist_ok=True)

            for path in mr_files:
                out_path = extraction_dir / (Path(path).stem + '.exc.tsv')
                try:
                    written = extract_excdata(path, out_path=out_path)
                    logger.info(f'wrote excitation data to {written}')
                except (ValueError, IndexError) as e:
                    logger.warning(f'{path}: no CASPT2/CAS-PDFT excitation data found ({e}); skipping')
                    continue

                ref_file = Path(opts.data_dir, "ref-excs.list")
                err_path = extract_excerr(written, ref_file)
                if err_path is not None:
                    logger.info(f'wrote {err_path}')

        results = {}
        for method in opts.methods:
            winners = run_method(method, parser, mr_calcs, ref_calcs, output_dir)
            results[method] = winners

            # Create symlinks from winners to excited-state data extraction files.
            if do_extract:
                create_excerr_symlinks(method, winners, mr_calcs, output_dir)

        if opts.generate_comparison_input:
            molcas_dir = output_dir / 'molcas_files'
            molcas_dir.mkdir(parents=True, exist_ok=True)
            written = generate_comparison_inputs(results, opts, molcas_dir)
            for path in written:
                logger.info(f'wrote timing input: {path}')

        all_results[keyword] = results

    return all_results
