### File description
In `src/`, you can find the following files.
- `GUGA_diag.py`: Evaluates the diagonal H matrix element of a CSF.
- `gen_spinfree_rdm.py`: Generates one- and two-body RDMs of a CSF.
- `IntegralClass.py`: Reads in an FCIDUMP file.
- `runmolcas.py`: Interfaces to run CSF-ROHF using OpenMolcas RASSCF module.
- **`settings.py`: To run `runmolcas.py`, you need to set your OpenMolcas build directory in this file.**

### Set OpenMolcas build directory
To use the interface, you need to set your OpenMolcas build directory in `settings.py`.

### Tutorial: SSG-CSF-ROHF Optimized Orbitals

This tutorial provides a step-by-step guide to obtaining SSG-CSF-ROHF optimized orbitals using the Fe dimer example from Reference [1].

Using the Python interface in this repository requires:
- [OpenMolcas](https://gitlab.com/Molcas/OpenMolcas) for orbital optimization  
- [NECI](https://github.com/fkfest/NECI_STABLE) for PT2-RDM generation  

> **Note:** The stochastic SplitGAS (SSG) method [2] is currently only available in an internal development version of NECI. Please contact us to request access.

All files used in this tutorial are located in `example/SSG-CSF-ROHF_tutorial`.

---
We begin with a set of localized and sorted high-spin ROHF orbitals `Fe2S2.SortOrb`.
The Fe 3d orbitals are site-separated: the first five orbitals belong to Fe 1 and
the last five orbitals belong to Fe 2.
Step 3 for PT2-RDMs generation is done with the SSG implementation [2] in NECI and
the other steps are donw with OpenMolccas.

#### Step 1: CSF-ROHF Orbital Optimization  
Directory: `1_CSF-ROHF/`

Run:
```bash
python ../../../src/runmolcas.py Fe2S2 1 1 1 1 1 2 2 2 2 2
```
`Fe2S2` is the input file name without extension.
The integers define the step-vector elements of the target CSF (|uuuuu ddddd>).
This command executes an OpenMolcas calculation via the Python interface and
it produces single-CSF-SCF variationally optimized orbitals.
For subsequent steps (2 and 4), we use the optimized orbitals from the final iteration before
canonicalization (`Fe2S2.IterOrb.25`).

#### Step 2: FCIDUMP Generation
Directory: `2_dumpgen/`

Run an OpenMolcas calculation using `Fe2S2.inp`.
The integral file is fed into NECI for PT2-RDM generation (step 3).

#### Step 3: PT2-RDM Generation
Directory: `3_RDMgen/`

Run NECI with the replica trick (`dneci`) using `input`.
`p1.list` defines the P space of the SSG calculation.

> **Note:** This step requires the internal development version of NECI with SSG support.

#### Step 4: SSG-CSF-ROHF Orbital Optimization
Directory: `4_SSG-CSF-ROHF/`

Run:
```bash
python ../../../src/runmolcas.py -u Fe2S2 1 1 1 1 1 2 2 2 2 2
```

The `-u` flag instructs the script to use RDM files from the working directory.
This produces SSG-CSF-ROHF orbitals and energy: `RASSCF energy for state  1 -5092.73230830`.
The energy is found in Table 3 of Reference [1].

> **Note:** The Python script supplies a fake energy to OpenMolcas to allow the iteration to proceed when the `-u` flag is used. This value appears in the output but should not be used for analysis.

### Other Example
Go to `example/` and `$ python ../src/runmolcas.py N2 1 1 1 2 2 2`.
`N2` is the Molcas input filenmae without extension, `1 1 1 2 2 2` is the CSF
you use for the ROHF optimization in the step-vector format.
This executes Molcas with `N2.inp` and feed RDMs and the RDM energy to Molcas
for every RASSCF iteration.

### Output
The python interface creates additional output files:
- `.pylog` logs the input parameters and the activites the interface does during the SCF iterations.
- `.iterdata` contains only the SCF iteration data of the calculation, as the OpenMolcas output contains additional information when
external RDMs are used (e.g., "echo $your_RDM_Energy ..." in between every iteration).

> **Note:** For SSG-CSF-ROHF calculations (Step 4 of the tutorial above), please only use `RASSCF_energy` (5th column) for the energy of the calculation.
`RDM_Energy` (last column) prints fake energies fed into OpenMolcas just to proceed iterations with the external RDM mode.

### References
- [1] Maru Song, Luca Bonfirraro, Ignacio Fdez. Galván, Roland Lindh, and Giovanni Li Manni,
"Spin-Adapted Restricted Open-Shell Hartree-Fock and Its Dynamic Correlation Extension",
Accepted for publication in J. Chem. Theory Comput. 2026,
https://doi.org/10.1021/acs.jctc.6c00379
- [2] Luca Bonfirraro, Oskar Weser, Maru Song, and Giovanni Li Manni,
"Stochastic-SplitGAS: A Quantum Monte Carlo Multi-Reference Perturbation Theory Based on the Imaginary-Time Evolution of Effective Hamiltonians",
J. Chem. Theory Comput. 2025, 21, 24, 12523–12544,
https://doi.org/10.1021/acs.jctc.5c01270
