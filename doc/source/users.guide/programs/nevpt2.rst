.. index::
   single: Program; NEVPT2
   single: NEVPT2

.. _UG\:sec\:nevpt2:

:program:`nevpt2`
=================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="NEVPT2">
            %%Description:
            <HELP>
            The NEVPT2 program computes the dynamic correlation correction to a DMRG-SCF energy
            according to the NEVPT2 theory developed originally by Angeli et al.
            Currently, NEVPT2 requires a DMRG reference wavefunction calculated with QCMaquis (with RASSCF or DMRGSCF module).
            </HELP>

NEVPT2 is a second-order perturbation theory with a CAS (or a CAS-like) reference wavefunction originally developed by Angeli et al. :cite:`Angeli_JChemPhys_Introduction_2001,Angeli_ChemPhysLett_Nelectron_2001,Angeli_JChemPhys_nelectron_2002,Angeli_JChemPhys_quasidegenerate_2004` In contrast to CASPT2, it uses a Dyall Hamiltonian :cite:`Dyall_JChemPhys_choice_1995` as the zeroth-order Hamiltonian and is therefore inherently free of intruder states and parameters such as the IPEA shift. NEVPT2 exists in two formulations -- the strongly- (SC-) and the partially-contracted NEVPT2 (PC-NEVPT2), which differ in the basis of the first-order wavefunction expansion.

The implementation in the :program:`NEVPT2` program is based on the original NEVPT2 implementation by Angeli et al. :cite:`Angeli_JChemPhys_nelectron_2002,Angeli_JChemPhys_quasidegenerate_2004`, with the implementation of the QCMaquis DMRG reference wave function and Cholesky decomposition for the two-electron integrals :cite:`Freitag_JChemTheoryComput_Multireference_2017`. For excited states both single-state and multi-state calculations with the QD-NEVPT2 approach :cite:`Angeli_JChemPhys_quasidegenerate_2004` are supported.

.. _UG\:sec\:nevpt2_dependencies:

Dependencies
------------
The :program:`NEVPT2` program needs the :file:`JOBIPH` file (or its HDF5 equivalent) with a reference wavefunction a from a :program:`RASSCF`/:program:`DMRGSCF` calculation. Currently, **only DMRG reference wavefunctions calculated with QCMaquis** are supported. Additionally, transformed MO integrals or Cholesky vectors from :program:`MOTRA` are required.

Optionally, four-particle reduced density matrices (and transition three-particle reduced density matrices for QD-NEVPT2 calculations) can be precalculated with QCMaquis in a massively parallel fashion and stored on disk. These QCMaquis calculations may be prepared and executed with the help of two scripts found in :file:`$MOLCAS/Tools/distributed-4rdm` directory, namely :file:`jobmanager.py` and :file:`prepare_rdm_template.sh`. The distributed RDM evaluation is strongly recommended for active spaces larger than 10-11 orbitals and is described in detail in :numref:`TUT:sec:nevpt2_distrdm`.

.. _UG\:sec\:nevpt2_input_files:

Input files
------------

.. class:: filelist

:file:`JobIph` or :file:`dmrgscf.h5`
  File containing information about the reference wavefunction.

:file:`ijkl.h5`
  Transformed integrals or Cholesky vectors, calculated by the :program:`MOTRA` program.

.. _UG\:sec\:nevpt2_output_files:

Output files
------------

.. class:: filelist

:file:`nevpt2.h5`
  File in HDF5 format, similar to RASSCF/DMRGSCF `dmrgscf.h5` files, containing the effective Hamiltonian for QD-NEVPT2 calculations (both strongly- and partially-contracted).

.. _UG\:sec\:nevpt2_input:

NEVPT2 input
------------

The :program:`NEVPT2` program is activated by ::

  &NEVPT2

The optional keywords supported by :program:`NEVPT2` are listed below.

.. class:: keywordlist

:kword:`STATES`
  Number of electronic states to calculate. Default: 1

  .. xmldoc:: <KEYWORD MODULE="NEVPT2" NAME="STATES" KIND="INT" LEVEL="BASIC">
              %%Keyword: STATES <basic>
              <HELP>
              Number of states to calculate. Default: 1
              </HELP>
              </KEYWORD>

:kword:`NOMS`
  Omit the QD-NEVPT2 calculation and perform single-state NEVPT2 calculations instead.

  .. xmldoc:: <KEYWORD MODULE="NEVPT2" NAME="NOMS" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: NOMS <basic>
              <HELP>
              Omit the (multi-state) QD-NEVPT2 calculation for multiple states.
              </HELP>
              </KEYWORD>

:kword:`MULT`
  Select specific states to perform QD-NEVPT2 calculation. Followed by a list of whitespace-separated state numbers, preceded by their total amount. Example: ``MULT=3 1 2 4`` for states 1, 2, 4 of a preceeding DMRG-SCF calculation of 4 roots (or more). ``MULT=ALL`` includes all states and is the default.

  .. xmldoc:: <KEYWORD MODULE="NEVPT2" NAME="MULT" KIND="CUSTOM" LEVEL="BASIC">
              %%Keyword: MULT <basic>
              <HELP>
              Select states for (multi-state) QD-NEVPT2 calculation.
              </HELP>
              </KEYWORD>

:kword:`FILE`
  Specify the path to a :file:`JobIph` or :file:`.h5` file with the reference wavefunction. By default, the reference wavefunction is read from :file:`JOBIPH`.

  .. xmldoc:: <KEYWORD MODULE="NEVPT2" NAME="FILE" KIND="STRING" LEVEL="BASIC">
              %%Keyword: FILE <basic>
              <HELP>
              Select JobIph or file with the reference wavefunction.
              </HELP>
              </KEYWORD>

:kword:`FROZEN`
  Specify the number of frozen orbitals. The number of frozen orbitals may be specified in two ways: if only one number :math:`n` is specified, then all orbitals from 1 to :math:`n` are frozen. Otherwise, it is possible to freeze specific orbitals with the :kword:`SELECT` keyword which follows the :kword:`FROZEN` keyword. In this case, the total number of frozen orbitals followed by the space-separated list of frozen orbitals must be entered. Note that if symmetry is used, the orbital numbering for all symmetries is still consecutive, e.g. the 1st orbital of symmetry 2 is has the number :math:`m+1` if there are :math:`m` orbitals in symmetry 1.

  If frozen orbitals are specified in :program:`MOTRA` input, they will be autodetected in :program:`NEVPT2` and there is no need to input them separately, so that this keyword is not needed.

  .. xmldoc:: <KEYWORD MODULE="NEVPT2" NAME="FROZEN" KIND="CUSTOM" LEVEL="BASIC">
              %%Keyword: FROZEN <basic>
              <HELP>
              Select frozen orbitals.
              </HELP>
              </KEYWORD>

:kword:`NOPC`
  Disable the PC-NEVPT2 calculation. If the option is not present (default), both SC-NEVPT2 and PC-NEVPT2 calculations are performed.

  .. xmldoc:: <KEYWORD MODULE="NEVPT2" NAME="NOPC" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: NOPC <basic>
              <HELP>
              Omit the partially-contracted NEVPT2 calculation.
              </HELP>
              </KEYWORD>

:kword:`SKIPK`
  Skip the calculation of Koopmans' matrices. Requires a file named :file:`nevpt.h5` obtained from a previous calculation in the scratch directory. May be useful to restart a previous crashed calculation if it crashed past the calculation of Koopmans' matrices, and may save some computational time, especially for large active spaces.

  .. xmldoc:: <KEYWORD MODULE="NEVPT2" NAME="SKIPK" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: SKIPK <basic>
              <HELP>
              Skip calculation of Koopmans' matrices and use them from a previous NEVPT2 calculation.
              </HELP>
              </KEYWORD>

:kword:`RDMRead`
  Do not calculate the 4-RDM, but rather read it from QCMaquis result files :file:`$Project.results_state.X.h5` for state ``X``. Useful if the previous calculation crashed but the 4-RDM evaluation step has succeeded. Do NOT use it if you are using the distributed 4-RDM calculation.

  .. xmldoc:: <KEYWORD MODULE="NEVPT2" NAME="RDMREAD" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: RDMRead <basic>
              <HELP>
              Read previously calculated 4-RDMs from QCMaquis result files instead of calculating it.
              </HELP>
              </KEYWORD>

:kword:`DISTributedRDM`
  Enable reading of the RDMs calculated with the distributed RDM evaluation script. This keyword should be followed by another line, which specifies the path to the folder with the calculation results. The 4-RDM will then be read from QCMaquis HDF5 files found in :file:`<path>/4rdm-scratch.<state>/parts/part-*/$Project.results_state.<state>.h5`. The distributed :math:`n`-RDM evaluation is described in the NEVPT2 program-based tutorial. If the tutorial is followed, the path should be :file:`$WorkDir`.

  .. xmldoc:: <KEYWORD MODULE="NEVPT2" NAME="DISTRIBUTEDRDM" KIND="STRING" LEVEL="BASIC">
              %%Keyword: DistributedRDM <basic>
              <HELP>
              Read RDMs calculated in a distributed fashion.
              </HELP>
              </KEYWORD>

.. _UG\:sec\:nevpt2_inputexample:

Input example
-------------

An input example for NEVPT2 may be found in :numref:`TUT:sec:nevpt2_run`.

.. xmldoc:: </MODULE>
