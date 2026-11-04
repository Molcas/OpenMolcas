# Dipole Moment Based Active Space (DM-AS) Chooser

Automated active-space selection for multireference calculations, based on
comparing candidate active space's dipole moments against reference dipole moments.

Details of the methods and their performance for various types of molecular systems can be found in:

> Kaufold, B. W.; Dong, S. S. Automated Active Space Selection with CASCI Dipole Moments, *J. Chem. Theory
> Comput.* 2026, DOI: [10.1021/acs.jctc.6c00473](https://doi.org/10.1021/acs.jctc.6c00473)

AND

> Kaufold, B. W.; Chintala, N.; Pandeya, P.; Dong, S. S. Automated Active Space Selection with Dipole Moments.
> *J. Chem. Theory Comput.* 2023, DOI: [10.1021/acs.jctc.2c01128](https://doi.org/10.1021/acs.jctc.2c01128).

## What it does

Given a set of candidate active-space calculations (from OpenMolcas CASCI/CASSCF, or a plain CSV) and a
set of reference dipole moments (from ORCA, or a plain CSV), this tool compares each candidate's
dipole moment(s) against each reference source, using one of six protocols, and reports which candidate
active space best reproduces the reference data.

## The eight protocols

| Protocol | Compares | Basis |
|---|---|---|
| `gdm-as` | Ground-state dipole magnitude | scalar |
| `edm-as` | Max absolute excited-state dipole difference, matched by state number | scalar |
| `d2dm-as` | Each side's independently-selected state of interest | scalar |
| `dir-gdm-as` | Ground-state dipole vector | vector |
| `dir-edm-as` | Max excited-state dipole vector-norm difference, matched by state number | vector |
| `dir-d2dm-as` | Each side's independently-selected state of interest, by vector norm | vector |
| `old-dir-gdm-as` | Ground-state dipole vector using 2023 formula | vector |
| `old-dir-edm-as` | Max excited-state dipole vector-norm difference, matched by state number, using 2023 formula | vector |

The difference between `dir-gdm-as` and `old-dir-gdm-as`, as well as between `dir-edm-as` and `old-dir-edm-as`, is the formula used to compare vector dipole moments between CASSCF/CASCI and the reference.  For `dir-gdm-as` and `dir-edm-as`, the "difference" between dipole moments is just the norm of the vector difference, as implemented in our most recent work.  For `old-dir-gdm-as` and `old-dir-edm-as`, the "difference" is a more involved calculation involving the projection of the CASSCF/CASCI dipole moment onto the reference dipole moment, as implemented in our 2023 work.  The formula used with `dir-gdm-as` and `dir-edm-as` is theoretically more robust but has only been tested with CASCI dipole moments, so we allow the user to use the old implementation if desired.

## Dependencies
 
- **`numpy`** : always required.
- [**`molextract`**](https://github.com/sdonglab/molextract) : only required if you're parsing real
  OpenMolcas/ORCA log files. You do not need to install it if you are providing CSV files for both candidates and references.
- **`pytest`** : only required to run the test suite; not needed to use the tool itself.

## Usage

```bash
python3 main.py \
  --methods gdm-as edm-as d2dm-as dir-gdm-as dir-edm-as dir-d2dm-as \
  --data-dir /path/to/data \
  -m casci casscf \
  -r b3lyp blyp camb3lyp lcwhpbe m06 m062x m06hf m06l pbe pbe0 r2scan0 r2scan50 r2scanh \
  --output-dir /path/to/results
```

### Arguments

- **`--methods`** : which protocol(s) to run (see table above).
- **`--data-dir`** : directory containing both candidate and reference files.
- **`--output-dir`** : where results are written (defaults to `./results`).
- **`-m` / `--mr-keywords`** : keyword(s) identifying candidate files by name. `-m casci casscf` finds
  every candidate file whose name contains `casci` or `casscf` as a delimited token.
- **`-r` / `--ref-keywords`** : keyword(s) identifying reference files by name. `-r b3lyp camb3lyp`
  builds two separate reference sources from whichever files contain each keyword as a delimited token.
- **`--debye-per-au`** : the conversion unit for converting Debye to atomic units. Default is 2.541765 (https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781118229101.app5)
- **`--hartree-to-ev`** : the conversion unit for converting Hartree to eV. Default is 27.211386245981 (https://physics.nist.gov/cgi-bin/cuu/Value?eqhrev)
- **`--extract-excdata`** : (optional) also run excitation-energy extraction/benchmarking (see below).
  Only applies to `-m` keywords that are CASSCF/CASPT2/CAS-PDFT calculations.
- **`--generate-comparison-input`** : (optional) also generate OpenMolcas timing-comparison input files for the
  winning active spaces (see below).

Each `-m` keyword is processed independently, with its own output subdirectory and its own set of
candidates. For example, `-m casci casscf` produces separate `results/casci/` and `results/casscf/` results.

## Candidate and reference file formats

Both candidates and reference sources can be provided either as real quantum-chemistry output files, or
as plain CSV files.

### OpenMolcas (.log) : candidates

A single log file containing a full CASCI/CASSCF calculation. Filenames must contain the calculation method 
(e.g. `casscf-mySystem.log`).

### ORCA (.out/.log/.orca) : reference sources

A single log file containing a ground-state (SCF) and/or an excited-state calculation. A
ground-state-only file has one dipole-moment block; an excited-state calculation has two.
The ground state calculation should be consistent throughout a single calculation method and system.
Filenames must contain the calculation method (e.g. `orca-pbe-mySystem.log`).

### CSV : either candidates or reference sources

```csv
state,dipole,x,y,z
0,1.9812,0.0,0.0,-1.9812
1,5.2057,0.15,-0.02,5.1980
```

- **`state`** is required for all protocols. `state=0` is the ground state; `1, 2, ...` are excited
  states.
- **`dipole`** is optional, only needed for the scalar protocols. 
- **`x`, `y`, `z`** are optional, only needed for the `dir-*` protocols.
- **All values (`dipole` and `x`/`y`/`z`) must be in Debye.**
- A single CSV can contain a candidate's or reference's full set of states (ground and excited together).
- Candidate CSV filenames must follow this naming format: `<keyword>-<electrons>-<orbitals>.csv`.
- Reference CSV filenames must include the keyword.

### Pooling multiple reference files for one source

If a reference source is built from several files (e.g. multiple TD-DFT calculations for the same
method), every file's states are pooled together. If more than one file provides a ground-state
(`state=0`) entry, all such entries must agree (on both magnitude and/or vector direction).

## Output

```
results/
└── casscf/                              # one subdirectory per -m keyword
    ├── active_space_selections.tsv      # winning active space per protocol and reference source
    ├── gdm-as.tsv                       # full comparison table for this protocol
    ├── edm-as.tsv
    ├── ...
    ├── BLYP_gdm                         # symlink: winning candidate's benchmark data (see below)
    ├── PBE_gdm
    └── excited_state_extractions/       # only present with --extract-excdata
        ├── casscf-06-08.exc.tsv
        └── casscf-06-08.excerr.tsv
    └── molcas_files/                    # only present with --generate-comparison-input
        ├── comparison-06-08.input
        └── comparison-10-14.input
```

Each protocol's `.tsv` file has one row per candidate active space, one column per reference source, and
each cell is the computed difference (empty if that comparison couldn't be computed, e.g. missing data).

`active_space_selections.tsv` records, for every protocol and reference source, which candidate active
space won (had the smallest difference).

## Generating OpenMolcas Input Files (optional)

The file `generate_active_space_input.py` is available at the root level to generate OpenMolcas input files for running active space calculations over a given range of active spaces.

The script can be run from the root with:
```bash
python3 generate_active_space_input.py \
  --output-dir path/to/desired/output/directory \   # defaults to cwd
  --coord your_coordinate_file.xyz \                # defaults to starting.xyz
  --basis "basis set in OpenMolcas format" \        # defaults to ANO-RCC-VTZP
  --nroots "# of roots to include in SA-CASCI" \    # defaults to 6
  --charge "molecular charge" \                     # defaults to 0
  --spin "the spin as (2S+1)" \                     # defaults to 1
  --nael-min "the minimum # of active space electrons to scan" \  # defaults to 6
  --nael-max "the maximum # of active space electrons to scan" \  # defaults to 14
  --nael-step "the step size between successive active space values" \  # defaults to 2
  --norb-max "the maximum # of active space orbitals to scan" \   # defaults to 14
  --norb-offset "the offset added to active-space-electrons/2 to get the minimum number of orbitals for a given active space" \  # defaults to 3
```

Running the output of this script will generate the OpenMolcas input files to use for the protocol.

## Excitation-energy extraction and benchmarking (optional)

Separate from active-space selection, `--extract-excdata` runs an additional benchmarking pipeline: it
extracts excitation energies (CASSCF/CASPT2/CAS-PDFT) from each candidate's OpenMolcas log, and, if a
`ref-excs.list` file (one reference excitation energy per line, no header) is present in `--data-dir`,
computes the absolute error against it for every candidate. This only works against real OpenMolcas
`.log` files, as CSV-sourced candidates have no excitation-energy data to extract.

For each protocol, the winning candidate's benchmarking data is symlinked at the top level of that
protocol's output directory (e.g. `BLYP_gdm -> excited_state_extractions/casscf-06-08.excerr.tsv`).

## Comparison input generation (optional)
 
`--generate-comparison-input` generates an OpenMolcas comparison input file (RASSCF/GRID_IT/CASPT2/MCPDFT/
RASSI) for every unique winning active space within each `-m` keyword group. If the same active space wins multiple times within one keyword group
(e.g. under both `gdm-as` and `edm-as`), only one input file is generated for it.
 
Generated files are named `comparison-<electrons>-<orbitals>.input` and written into that keyword's own
output subdirectory, alongside its other results.  

### Additional arguments (only relevant with `--generate-input`)
 
- **`--coord`** : path to the system's Cartesian coordinates file (default: `starting.xyz`).
- **`--basis`** : basis set in OpenMolcas format (default: `ANO-RCC-VTZP`).
- **`--integrals`** : path to pre-computed two-electron integrals (default: `$CurrDir/2e-integrals/2e-integrals.OrdInt`).
- **`--fileorb`** : path to starting orbitals for RASSCF (default: `$CurrDir/uhf/uhf.UnaOrb`).
- **`--nroots`** : number of roots for SA-CASSCF, also sets the GRID_IT loop ranges (default: `6`).
- **`--charge`** : molecular charge (default: `0`).
- **`--spin`** : spin as (2S+1) (default: `1`).
- **`--cimx`** : maximum RASSCF CI iterations (default: `200`).
- **`--caspt2-maxiter`** : maximum CASPT2 iterations (default: `10000`).
- **`--ksdft`** : MC-PDFT functional (default: `T:PBE`).

## Running the tests

```bash
python3 -m pytest tests/
```
Note that running the tests requires [**`molextract`**](https://github.com/sdonglab/molextract) to be installed.

## Project structure

```
main.py                             # CLI entry point (argument parsing only)
generate_active_space_input.py      # To generate OpenMolcas .input files for running active space calculations
select_active_spaces.sh             # How to call the protcol if submitting to an HPC cluster
src/
  workflow_orchestrator.py          # workflow: file discovery, candidate/reference resolution, per-protocol runs
  calc_base.py                      # MultiRefCalc, SingleRefCalc
  parsers.py                        # parse_orca_log, parse_csv_states, parse_mr_log
  extract_excited.py                # extract_excdata, extract_excerr (benchmarking pipeline)
  choose.py                         # picks the winning active space for each reference source
  utils.py                          # shared structures and constants (ElectronicState, Vector)
  generate_molcas_files.py          # to generate comparison .input files for the selected active spaces
  as_methods/
    base.py                                   # Table (abstract base): filter_input, filter_ground, filter_excited,
                                              #     check_consistent_ground_states, build, write_tsv
    gdm.py, edm.py, d2dm.py                   # scalar protocols
    dir_gdm.py, dir_edm.py, dir_d2dm.py       # vector protocols
    old_dir_gdm.py, old_dir_edm.py            # vector protocols (2023 version)
tests/
  orchestrator_test.py
  extract_excited_test.py
  csv_test.py
  choose_test.py
  generate_molcas_files_test.py
  method_tests/
    gdm_test.py, edm_test.py, d2dm_test.py
    dir_gdm_test.py, dir_edm_test.py, dir_d2dm_test.py
    old_dir_gdm_test.py, old_dir_edm_test.py
```
## Authors
Copyright (C) 2026, Meredith E. Wanta, Benjamin W. Kaufold, and Sijia S. Dong, Northeastern University

