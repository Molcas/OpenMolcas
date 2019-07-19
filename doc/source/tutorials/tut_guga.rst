.. index::
   single: Program; GUGA
   single: GUGA

.. _TUT\:sec\:guga:

:program:`GUGA` --- A Configuration Interaction Coupling Coefficients Program
=============================================================================

Several of the Configuration Interaction (CI) modules in |molcas| use
the :program:`guga` module to compute the CI coupling coefficients.
We continue our water calculations using the input file shown in
the input below. The :kword:`TITLe` keyword behaves
in a similar fashion as described in previous modules.
There are several compulsory keywords of the :program:`guga` module. The
number of electrons to be correlated is specified using the
:kword:`ELECtrons` keyword. We are correlating the valence electrons.
The spin state is specified using the :kword:`SPIN` keyword.

.. index::
   single: GUGA; Input

Sample input requesting the the GUGA module to calculate the coupling
coefficients for neutral triplet water in :math:`C_{2v}` symmetry with six electrons
in the active space: ::

  &GUGA
  Title= GUGA for C2v Water
  Electrons= 8; Spin= 3
  Inactive= 1 0 0 0; Active= 2 2 0 1
  CIAll= 1

The keywords :kword:`CIALl` and :kword:`REFErence` are mutually
exclusive. We specify :kword:`CIALl` which will calculate the
energy using all possible references functions that can be constructed
using the input set of occupation numbers of the active orbitals regardless of
the spin coupling (all configurations used to build the corresponding CASSCF
wave function). Specific selected references can be chosen using the
:kword:`REFErence` keyword. Either the :kword:`ACTIve` or :kword:`INACtive`
keyword should be used for a meaningful calculation. The default for both
keywords is zero for all symmetries. These keywords function in a similar
fashion to these in the :program:`RASSCF` program module. The :kword:`INACtive`
keyword specifies the orbitals that are fully occupied in each symmetry
in all the reference functions and the :kword:`ACTIve` keyword
specifies the orbitals that may have varying occupations in all references.
The selection of :kword:`INACtive` orbitals in the input above
is forcing the bonding sp hybrid orbital to remain fully occupied in all
reference states.

:program:`GUGA` Output
----------------------

The :program:`GUGA` section of the output lists the possible configurations
in the active space. There are nine possible triplet configurations of
six electrons in five orbitals. Apart from the various types of orbital in each
symmetry the :program:`GUGA` section of the output also gives the number of
states that will coupled with various states. There are no input files for the
:program:`GUGA` module but the calculated coupling coefficients are stored in
:file:`CIGUGA`.
