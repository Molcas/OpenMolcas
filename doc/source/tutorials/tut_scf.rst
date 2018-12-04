.. index::
   single: Program; SCF
   single: SCF

.. _TUT\:sec\:scf:

:program:`SCF` --- A Self-Consistent Field program and Kohn--Sham DFT
=====================================================================

The simplest *ab initio* calculations possible use the Hartree--Fock
(HF) Self-Consistent Field (SCF) method with the program name :program:`SCF` in
the |molcas| suite. It is possible to calculate the HF energy once we have
calculated the integrals using the :program:`SEWARD` module, although |molcas|
can perform a direct SCF calculation in which the two-electron integrals are
not stored on disk. The |molcas| implementation performs a closed-shell (all
electrons are paired in orbitals) and open-shell (Unrestricted Hartree--Fock)
calculation. It is not possible to perform an Restricted Open-shell Hartree--Fock (ROHF)
calculation with the :program:`SCF`. This is instead done using the program
:program:`RASSCF`. The :program:`SCF` program can also be used to perform
calculations using Kohn Sham Density Functional Theory (DFT).

The :program:`SCF` input for a Hartree--Fock calculation of a water
molecule is given in :numref:`block:scf_input`
which continues our calculations on the water molecule.

There are no compulsory keywords following the program name, ``&SCF``. If no input
is given the program will compute the SCF energy for a neutral molecule with the
orbital occupations giving the lowest energy. Here, we have used the following
input: the first is :kword:`TITLe`. As
with the :program:`SEWARD` program, the first line following the keyword is
printed in the output.

.. index::
   single: SCF; Input
   single: SCF; Occupied

No other keyword is required for a closed-shell calculation. The program
will find the lowest-energy electron configuration compatible with the
symmetry of the system and will distribute the orbitals accordingly.
In complex cases the procedure may fail and produce a higher-lying configuration.
It is possible to use the keyword :kword:`OCCUpied`
which specifies the number of occupied orbitals in each symmetry grouping
listed in the :program:`GATEWAY` output and given in
:numref:`block:Irreducible`, forcing the method to converge to the specified
configuration. The basis label and type give an
impression of the possible molecular orbitals that will be obtained in
the SCF calculation. For example, the first basis function in the :math:`a_1`
irreducible representation is an s type on the oxygen indicating the
oxygen 1s orbital. Note, also, that the fourth basis function is centered on
the hydrogens, has an s type and is symmetric on both hydrogens as
indicated by both hydrogens having a phase of 1, unlike the sixth basis function
which has a phase of 1 on center 2 (input H1) and |-|\ 1 on center 3
(generated H1).
As an alternative you can use the keyword :kword:`Charge` with parameters 0 and
1 to indicate a neutral molecule and optimization procedure 1 that searches for
the optimal occupation.

.. code-block:: none
   :caption: Sample input requesting the :program:`SCF` module to calculate the ground Hartree--Fock energy for a neutral water molecule in :math:`C_{2v}` symmetry.
   :name: block:scf_input

   &SCF
   Title= Water - A Tutorial. The SCF energy of water is calculated using C2v symmetry
   End of Input

.. index::
   single: Symmetry; Adapted basis functions

.. code-block:: none
   :caption: Symmetry adapted Basis Functions from a :program:`GATEWAY` output.
   :name: block:Irreducible

             Irreducible representation : a1
             Basis function(s) of irrep: z

   Basis Label        Type   Center Phase
     1   O1           1s        1     1
     2   O1           2s        1     1
     3   O1           2p0       1     1
     4   H1           1s        2     1      3     1

             Irreducible representation : b1
             Basis function(s) of irrep: x, xz, Ry

   Basis Label        Type   Center Phase
     5   O1           2p1+      1     1
     6   H1           1s        2     1      3    -1

             Irreducible representation : b2
             Basis function(s) of irrep: y, yz, Rx

   Basis Label        Type   Center Phase
     7   O1           2p1-      1     1

.. Note: This includes a nbsp character

We have ten electrons to ascribe to five orbitals to describe a
neutral water molecule in the ground state. Several
techniques exist for correct allocation of electrons. As a test of
the electron allocation, the energy obtained should be the same with
and without symmetry.
Water is a simple case, more so when
using the minimal basis set. In this case, the fourth irreducible
representation is not listed in the :program:`GATEWAY` output as there
are no basis functions in that representation.

.. index::
   single: SCF; Open-shell cases – Unrestricted Kohn–Sham DFT

To do a UHF calculation, the keyword :kword:`UHF` must be specified.
To force a specific occupation for alpha and beta orbitals
In this keyword :kword:`OCCNumbers` has to be used with two entries, one
for alpha and beta occupied orbital. It is possible to use UHF
together with keyword :kword:`Charge` or :kword:`Aufbau`, in this case
you have to specify a keyword :kword:`ZSPIN` set to the
difference between alpha and beta electrons.

If you want to do an UHF calculation for a closed shell system, for example,
diatomic molecule with large interatomic distance, you have to specify keyword
:kword:`SCRAMBLE`.

To do the Density Functional Theory calculations, keyword :kword:`KSDFT` followed
in the next line by the name of the available functional as listed in the input
section is compulsory. Presently following Functional Keywords are available:
BLYP, B3LYP, B3LYP5, HFB, HFS, LDA, LDA5, LSDA, LSDA5, SVWN, SVWN5, TLYP, PBE, PBE0,
M06, M06HF, M062X, M06L.
The description of functional keywords and the functionals is defined in the
section DFT Calculations :ref:`UG:sec:dft`.

The input for KSDFT is given as, ::

  KSDFT= B3LYP5

In the above example B3LYP5 functional will be used in KSDFT calculations.

Running :program:`SCF`
----------------------

Performing the Hartree--Fock calculation introduces some important
aspects of the transfer of data between the |molcas| program modules.
The :program:`SCF` module uses the integral files computed by
:program:`SEWARD`. It produces a orbital file with the symbolic name
:file:`SCFORB` which contains all the MO information. This is then
available for use in subsequent |molcas| modules. The
:program:`SCF` module also adds information to the :file:`RUNFILE`.
Recall that the :program:`SEWARD` module produces two integral files
symbolically linked to :file:`ONEINT` and :file:`ORDINT` and actually
called, in our case, :file:`water.OneInt` and :file:`water.OrdInt`,
respectively (this is for non-Cholesky-type calculations only).
Because the two integral files are present in
the working directory when the :program:`SCF` module is performed, |molcas|
automatically links them to the symbolic names.

If the integral files were not deleted in a previous calculation
the :program:`SEWARD` calculation need not be repeated. Furthermore,
integral files need not be in the working directory if they are linked
by the user to their respective symbolic names. Integral files,
however, are often very large making it desirable to remove them after the
calculation is complete. The linking of files to their symbolic names
is useful in other case, such as input orbitals.

.. index::
   single: SCF; LumOrb
   single: SCF; Input orbitals
   single: SCF; Convergence problems

If nothing else is stated, the :program:`SCF` program will use the guess orbitals
produced by :program:`SEWARD` as input orbitals with the internal name
:file:`GUESSORB`. If one wants to use any other input orbitals for the
:program:`SCF` program the option :kword:`LUMOrb` must be used. The
corresponding file should be copied to the internal file :file:`INPORB`. This
could for example be an orbital file generated by an earlier SCF calculation,
:file:`$Project.ScfOrb`. Just copy or link the file as :file:`INPORB`.

.. index::
   single: SCF; Output

:program:`SCF` Output
---------------------

The :program:`SCF` output includes the title from the input as well as
the title from the :program:`GATEWAY` input because we used the integrals
generated by :program:`SEWARD`. The output also contains the cartesian
coordinates of the molecule and orbital specifications including the
number of frozen, occupied and virtual (secondary) orbitals in each
symmetry. This is followed by details regarding the :program:`SCF`
algorithm including convergence criteria and iteration limits. The
energy convergence information includes the one-electron, two-electron,
and total energies for each iteration. This is followed by the final
results including the final energy and molecular orbitals for each
symmetry.

The Density Functional Theory Program gives in addition to the above,
details of grids used, convergence criteria, and name of the functional used.
This is followed by integrated DFT energy which is the functional contribution
to the total energy and the total energy including the correlation.
This is followed results including the Kohn Sham orbitals for each symmetry.

.. index::
   single: SCF; Orbitals
   single: SCF; Orbital energies
   single: Orbital energies

The molecular orbital (MO) information lists the orbital energy, the
electron occupation and the coefficients of the basis functions
contributing to that MO. For a minimal basis set, the basis functions
correspond directly to the atomic orbitals. Using larger basis sets
means that a combination of the basis functions will be used for each
atomic orbital and more so for the MOs.
The MOs from the first symmetry species are
given in :numref:`block:water_MOs`. The first MO has an energy of
-20.5611 hartree and an occupation of 2.0. The major
contribution is from the first basis function label "``O1  1s``"
meaning an s type function centered on the oxygen atom. The
orbital
energy and the coefficient indicates that it is the MO based largely
on the oxygen 1s atomic orbital.

.. code-block:: none
   :caption: Molecular orbitals from the first symmetry species of a calculation of water using :math:`C_{2v}` symmetry and a minimal basis set.
   :name: block:water_MOs

   Molecular orbitals for symmetry species 1: a1

      Orbital        1         2         3         4
      Energy    -20.5611   -1.3467    -.5957     .0000
      Occ. No.    2.0000    2.0000    2.0000     .0000

    1 O1  1s      1.0000    -.0131    -.0264    -.0797
    2 O1  2s       .0011     .8608    -.4646    -.7760
    3 O1  2p0      .0017     .1392     .7809    -.7749
    4 H1  1s      -.0009     .2330     .4849    1.5386

The second MO has a major contribution from the second oxygen ``1s``
basis function indicating a mostly oxygen ``2s`` construction.
Note that it is the
absolute value of the coefficient that determines it importance. The
sign is important for determining the orthogonality of its orbitals and
whether the atomic orbitals contributions with overlap constructively
(bonding) or destructively (anti-bonding).
The former occurs in this MO as indicated by the positive sign on the
oxygen ``2s`` and the hydrogen ``1s`` orbitals, showing a bonding
interaction between them.
The latter occurs in the third MO, where the relative sign is reversed.

The third MO has an energy of
-0.5957 hartree and major contributions from the second oxygen
``1s`` basis function, the oxygen ``2p0`` basis function and the
hydrogen ``1s`` basis functions which are symmetrically situated on
each hydrogen (see :numref:`block:Irreducible`). The mixing of the
oxygen ``2s`` and ``2p0`` basis functions leads to a hybrid
orbital that points away from the two hydrogens, to which it is
weakly antibonding.

A similar analysis of the fourth orbital reveals that it is the
strongly anti-bonding orbital partner to the third MO. The oxygen ``2p0``
basis function is negative which reverses the overlap characteristics.

The molecular orbital information is followed by a Mulliken charge
analysis by input center and basis function. This provides a measure
of the electronic charge of each atomic center.

Towards the end of the :program:`SCF` section of the |molcas| output
various properties of the molecule are displayed. By default the
first (dipole) and second cartesian moments and the quadrupoles are displayed.

.. The
   inclusion of the :kword:`FLDG` keyword (with zero (0))
   with cause the electric field gradients at each atomic center to be calculated
   and displayed.
   There are several other properties that can be calculated
   in this fashion using the variational |molcas| programs --- :program:`SCF`
   and :program:`RASSCF` when producing a CASSCF wave function.

:program:`SCF` --- Basic and Most Common Keywords
-------------------------------------------------

.. class:: keywordlist

:kword:`UHF`
  Unrestricted Hartee Fock or unrestricted DFT calculation

:kword:`KSDFt`
  DFT calculations, with options: BLYP, B3LYP, B3LYP5, HFB, HFS,
  LDA, LDA5, LSDA, LSDA5, SVWN, SVWN5, TLYP, PBE, PBE0, M06, M06HF, M062X, M06L

:kword:`CHARge`
  Net charge of the system (default zero)

:kword:`ZSPIn`
  Difference between :math:`\alpha` and :math:`\beta` electrons

:kword:`Occupied`
  Specify the orbital occupations per irreps
