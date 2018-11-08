.. index::
   single: Program; MPSSI
   single: MPSSI

.. _UG\:sec\:mpssi:

:program:`mpssi`
================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="MPSSI">
            %%Description:
            <HELP>
            %%Description:
            In analogy to the RASSI program, the MPSSI program calculates overlaps, and matrix
            elements of one-electron operators, and of the electronic Hamiltonian,
            over a basis of matrix-product state (MPS) wave functions calculated by QCMaquis,
            which may each have its own
            independent set of orbitals. Energies and matrix elements are
            computed also for the non-interacting linear combinations of states,
            i.e., doing a limited CI using the MPS states as a non-orthogonal basis.
            MPSSI can be used to compute dipole oscillator strengths,
            spin-orbit interaction matrix elements as well as, for example, transition dipole
            moments for spin-orbit eigenstates.
            </HELP>

The
:program:`MPSSI` (MPS :index:`State Interaction`) program forms overlaps and
other :index:`matrix
elements <single: Matrix elements; DMRGSCF>` of the Hamiltonian and other operators
over a wave function basis, which consists of matrix-product state (MPS) wave functions,
each with an individual set of orbitals.
Following the philosophy of the :program:`RASSI` program, it is a generalized state-interaction approach
for both nonorthogonal and orthonormal spinfree
MPS wave functions which enables the evaluation of arbitrary one- and two-particle
transition matrix elements as well as, for example, matrix elements of the
spin-orbit coupling operator. For instance, diagonalization of the spin-orbit Hamiltonian
matrix yields spin-orbit coupled wave functions as linear combinations of the
uncoupled, spin-pure MPS states. The latter can (but do not have to) be obtained as
results from one or several DMRG-SCF orbital optimization calculations (see :program:`DMRGSCF`).

Following the work of Malmqvist :cite:`Malmqvist:89`, the central element of the MPS-SI algorithm
is the transformation of the bra and ket MPS wave functions to a biorthonormal basis
representation. It is important to note that the latter transformation is not needed if
the MPS wave functions considered for state interaction share a common MO basis.
In this particular case, the MPS-SI program directly proceeds with the calculation of
the reduced (transition) one- and two-particle density matrices. We emphasize that
our approach is applicable to the general case with MPS wave functions built from
mutually nonorthogonal molecular orbital bases. It therefore provides the desired
flexibility to find the best individual molecular orbital basis to represent wave functions
of different spin and/or spatial symmetry. After solving a generalized eigenvalue
equation of the form

.. math:: Hc = ESc
   :label: egv

with the Hamiltonian matrix :math:`H` expressed in the basis of the DMRG-SCF MPS wave
functions and the wave function overlap matrix :math:`S`, a set of fully orthogonal and noninteracting
states are obtained as linear combinations of the DMRG-SCF MPS wave
functions with the expansion coefficients given by :math:`c` in Eq. :eq:`egv`.

Apart from computing oscillator strengths, overlaps and Hamiltonian
matrix elements can be used to compute :index:`electron transfer rates <single: Electron transfer rate>`, or
to form :index:`quasi-diabatic states <single: Quasi-diabatic states>` and reexpress matrix elements over a
basis of such states.

Moreover, it is possible to “dress" the diagonal elements of the Hamiltonian in
Eq.  :eq:`egv` for MPS-SI by adding a correlation-correction term obtained, for example,
from a preceding NEVPT2 calculation (see Section 6), by either using the HDIAG keyword
within the RASSI module or provide the nevpt2.h5 wave function file as input


.. index::
   pair: Dependencies; MPSSI

.. _UG\:sec\:mpssi_dependencies:

Dependencies
------------

The :program:`MPSSI` program needs one or more :file:`dmrgscf.h5` files produced
by the :program:`DMRGSCF` program (or if :program:`MPSSI` is running subseuqently after a :program:`NEVPT2` calculation
one or more :file:`nevpt2.h5` files). Also, it needs a :file:`ONEINT` file from
:program:`SEWARD`, with overlap integrals and any one-electron
property integrals for the requested matrix elements. If Hamiltonian
matrix elements are used, also the :file:`ORDINT` file is needed.

For further information see the description of the :program:`RASSI`.

.. index::
   pair: Files; MPSSI

.. _UG\:sec\:mpssi_files:

Files
-----

Input files
...........

.. class:: filelist

:file:`ORDINT*`
  Ordered two-electron integral file produced by the :program:`SEWARD`
  program. In reality, this is up to 10 files in a multi-file system,
  named ORDINT, ORDINT1,...,ORDINT9. This is necessary on some platforms
  in order to store large amounts of data.

:file:`ONEINT`
  The one-electron integral file from :program:`SEWARD`

:file:`dmrgscf.h5`
  A number of :file:`dmrgscf.h5` files from different :program:`DMRGSCF` jobs.

Output files
............

.. class:: filelist

:file:`SIORBnn`
  A number of files containing natural orbitals, (numbered sequentially as
  :file:`SIORB01`, :file:`SIORB02`, etc.)

:file:`BRAORBnnmm`, :file:`KETORBnnmm`
  A number of files containing binatural orbitals for the transition between
  states nn and mm.

:file:`TOFILE`
  This output is only created if :kword:`TOFIle` is given in the input.
  It will contain the transition density matrix computed by :program:`MPSSI`.
  Currently, this file is only used as input to :program:`QmStat` (NOT TESTED!).

:file:`EIGV`
  Like :file:`TOFILE` this file is only created if :kword:`TOFIle` is given
  in the input. It contains auxiliary information that is picked up
  by :program:`QmStat` (NOT TESTED!).

  .. :file:`UNSYM`
       The derivative of the transition dipole moment desymmetrized.

.. index::
   pair: Input; MPSSI

.. _UG\:sec\:mpssi_input:

Input
-----

This section describes the input to the
:program:`MPSSI` program in the |molcas| program system,
with the program name: ::

  &MPSSI

When a keyword is followed by additional mandatory lines of input,
this sequence cannot be interrupted by a comment line. The first 4
characters of keywords are decoded. An unidentified keyword makes the
program stop. Note that :program:`MPSSI` shares **ALL** keywords with :program:`RASSI` which do **NOT** request CI-type
quantities. Below is just a list of additional keywords available for enabling the effective Hamiltonian from a
preceeding :program:`NEVPT2` calculation, in order to achieve a state-dressing.

.. index::
   pair: Keywords; MPSSI

Keywords
........

.. class:: keywordlist

:kword:`QDSC`
  Enable the effective Hamiltonian from a quasi-degenerate (QD) multi-state strongly-contracted i(SC) :program:`NEVPT2`
  calculation.

  .. xmldoc:: <GROUP MODULE="MPSSI" NAME="QDSC" APPEAR="effective Hamiltonian" KIND="BLOCK" LEVEL="ADVANCED">
              %%Keyword: QDSC <advanced>
              <HELP>
              Enable the QDSC effective Hamiltonian in MPSSI.
              </HELP>

:kword:`QDPC`
  Enable the effective Hamiltonian from a quasi-degenerate (QD) multi-state partially-contracted (PC) :program:`NEVPT2`
  calculation.

  .. xmldoc:: <GROUP MODULE="MPSSI" NAME="QDPC" APPEAR="effective Hamiltonian" KIND="BLOCK" LEVEL="ADVANCED">
              %%Keyword: QDPC <advanced>
              <HELP>
              Enable the QDPC effective Hamiltonian in MPSSI.
              </HELP>

Input example
.............

::

  &MPSSI &END
  NrofJobIphs
  1 2           --- 1 JobIph (actually an .h5 file) - 2 states to be read
  1 2           --- which roots from the .h5 file.
  FILE
  1
  n2+.dmrgscf.h5
  omega
  SPIN
  EPRG
  1.0
  MEES
  PROPerties
    3
  'AngMom' 1
  'AngMom' 2
  'AngMom' 3
  * This input will compute spinfree and spin-orbit igenstates in the space
  * spanned by the 2 input functions

.. xmldoc:: </MODULE>
