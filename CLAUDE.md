# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

OpenMolcas is a quantum chemistry package (mostly modern Fortran, `.F90` with preprocessing, plus some C and Python), focused on multiconfigurational methods (CASSCF, CASPT2, ...). LGPL 2.1.

## Build

Out-of-source CMake build (the build directory doubles as the runtime `$MOLCAS` directory):

```
git submodule update --init External/lapack   # only if not using MKL/OpenBLAS
mkdir build && cd build
cmake ../OpenMolcas            # options: -DMPI=ON -DOPENMP=ON -DHDF5=ON -DBOUNDS=ON -DBIGOT=ON -DCMAKE_BUILD_TYPE=Debug ...
make -j
./pymolcas -setup              # first install
```

`export MOLCAS=<build dir>` to run calculations or the helper scripts below. Useful options: `BOUNDS` (gfortran bounds checking), `BIGOT` (warnings as errors, as in CI), `GCOV`, `BUILD_TESTING` (unit tests, ON by default). Optional external interfaces (DMRG, NECI, WFA, MSYM, ...) are off by default and usually need a submodule in `External/`.

## Tests

Two kinds:

- **Verification tests** (the main regression suite): `test/<group>/NNN.input` (groups: `standard`, `additional`, `exhaustive`, `grayzone`, ...). Each file is a Molcas input followed by an embedded `>>FILE checkfile ... >>EOF` block holding reference values (e.g. `E_RASSCF="..."/8` = value/tolerance digits). Run from the build directory:
  ```
  ./pymolcas verify                     # default set; grayzone failures are expected
  ./pymolcas verify standard:000        # single test (or just: 000)
  ./pymolcas verify standard:050-070
  ./pymolcas verify --module rasscf     # tests that use a module
  ./pymolcas verify --grep WORD         # tests containing a keyword
  ./pymolcas verify -k -d standard:000  # keep work dir, print output to terminal
  ./pymolcas verify --failed            # rerun last failures
  ./pymolcas verify --generate <test>   # regenerate checkfile (only when a change in results is intended)
  ```
  New test files must carry a checkfile whose `* Molcas version <commit>` is already in the branch history (enforced by `sbin/check_tests` in the pre-commit hook): commit the code first, rebuild, then `--generate`. A test without a checkfile stops at the first module unless run with `--pass`.
- **CSF/SD verification** (CSF/SD implementation, tested through the `ci_interfaces` routines): `test/csf_sd/` (light, seconds; run after any change in `src/sguga_util/`; tests 0xx are without symmetry, 1xx with symmetry) and `test/csf_sd_heavy/` (minutes). Both carry a `.csf_sd` dot file and no default-group dot files, so they never run by default:
  ```
  ./pymolcas verify csf_sd
  MOLCAS_TIMELIM=3600 ./pymolcas verify -j 2 csf_sd_heavy  # verify kills tests after 600 s by default
  ./pymolcas verify .csf_sd                                # both
  ```
  The inputs set `>export MOLCAS_CI_VERIFY=...`, which makes RASSCF call `SG_Test_EpqErs_Ref` and abort on any mismatch against the reference codes: `YES` checks `<I|Epq|J>` and `<I|EpqErs|J>`, `EPQ` only `<I|Epq|J>` (used where `MatEl_EpqErs` is not yet debugged for symmetry, e.g. the H8 cube tests `csf_sd:010,109-115`). Without the variable, ordinary RASSCF runs skip the checks. The checks scale steeply with the number of CSFs (a few hundred CSFs already take minutes), so keep new tests small or put them in `csf_sd_heavy`.
  `test/csf_sd_diag/` (no `.csf_sd` file) holds diagnostic tests that are **expected to fail** until a known bug is fixed (e.g. `csf_sd_diag:101`, the H4 tetrahedron in D2, `<I|Epq|J>` mismatch for upper levels in different irreps). Their checkfiles were generated with `MOLCAS_CI_VERIFY=NO` (the energies do not depend on the checks) and the `>export` line then switched back; once one passes, move it to `test/csf_sd`. Note that the RASSCF sigma vector comes from Lucia/Faroald, not from SGUGA, so correct energies do not validate `SG_Epq_Psi`; the SGUGA sigma is only compared against Lucia when compiled with `_SGUGA_VERIFY_` (`ci_interfaces.F90`).
- **Unit tests** (FRUIT-based, `unit_tests/`): built with the project, run via `ctest` in the build directory (`ctest -R <name>` for one).

## Style checks

`sbin/check_style` (with `MOLCAS` pointing at the source tree, or via the pre-commit hook installed by `sbin/install_hooks.sh`) enforces, for Fortran sources: ASCII only, no tabs, no trailing blanks, no `write(*,...)`, no `double precision`, no `STOP`, `PRINT`, raw `OPEN`, `MATMUL`, `FORALL`, `character*n`, numeric `GOTO`, `ERR=`/`END=`, no bare BLAS/LAPACK calls (use the wrappers), no `_DEBUG_` flag. The hook also runs `sbin/copyright` and `sbin/check_tests`, and blocks commits to `master`. Formatting follows `Tools/fprettify/fprettifyrc` (2-space indent, 132-column lines, lowercase keywords).

## Fortran conventions

- Kinds and units come from `Definitions` (`use Definitions, only: wp, iwp, u6`); write `real(kind=wp)`, `integer(kind=iwp)`, output to `u6`.
- Numeric constants from `Constants` (`Zero`, `One`, `Two`, `Half`, ...).
- Memory via `stdalloc` (`mma_allocate`/`mma_deallocate`, `mma_maxDBLE`), not plain `allocate`.
- Explicit `use ..., only:` lists; `implicit none`. Debug code is guarded by `#define _DEBUGPRINT_` / `#ifdef _DEBUGPRINT_` (not `_DEBUG_`).
- Every new file needs the OpenMolcas LGPL header with copyright lines (one author per line, all years); add new authors to `CONTRIBUTORS.md` (each line ends with two spaces, sorted by last name). See `CONTRIBUTING.md`.

## Architecture

- `src/<name>/` without `_util` suffix is a **program module** (e.g. `seward`, `scf`, `rasscf`, `caspt2`, `mclr`, `alaska`, `slapaf`), invoked as `&NAME` in an input file. `src/<name>_util/` directories are **utility libraries**. The top-level `CMakeLists.txt` discovers both by globbing `src/` (so a new directory is picked up automatically) and links all utils into `libmolcas`; per-program extra deps are `<prog>_deplibs` there.
- Each directory's `CMakeLists.txt` lists `sources` explicitly (new files must be added) and includes `cmake/prog_template.cmake` or `cmake/util_template.cmake`. In a util, files defining modules used by other utils must be listed in `modfile_list`.
- `src/Include/` holds shared `.fh` include files (e.g. `intent.fh`, interface headers); `src/Driver/` holds the `*.prgm.src` files declaring each program's input/output files for the driver; the driver itself is `Tools/pymolcas`.
- Programs communicate through files in the work directory, chiefly the **RunFile** (`runfile_util`: `Get_*`/`Put_*` routines) plus integral/orbital files (`OneInt`, `OrdInt`, `JobIph`, HDF5 `.h5`). A typical chain is `gateway` → `seward` (integrals, Cholesky/RI) → `scf` → `rasscf` → `caspt2`/`rassi`/`mclr`/`alaska` → `slapaf`.
- CI solvers for RASSCF live in several utils (`lucia_util`, `molcas_ci_util`, `gugaci`/`guga_util`, and on this branch `sguga_util`, a symmetric-group GUGA implementation with the `SGUGA` module in `sguga.F90`).
- `basis_library/` and `data/` are runtime data; `doc/` is the Sphinx manual; `Tools/` holds optional tools (built with `-DTOOLS=ON`); `sbin/` holds dev scripts (`help_func`, `find_sources`, etc., also reachable as `pymolcas <script>`).

## Commits

Short title line, blank line, explanation (see `CONTRIBUTING.md`). Upstream is GitLab (`gitlab.com/Molcas/OpenMolcas`), contributions go via merge requests.
