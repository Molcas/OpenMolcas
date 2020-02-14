.. index::
   single: Program; GENANO
   single: GENANO

.. _UG\:sec\:genano:

:program:`genano`
=================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="GENANO">
            %%Description:
            <HELP>
            This program is used to construct ANO type basis sets.
            </HELP>

:program:`GENANO` is a program for
determining the contraction coefficients for
generally contracted basis sets :cite:`Raffenetti:73`.
They are determined by diagonalizing a density matrix,
using the eigenvectors (natural orbitals) as
the contraction coefficients, resulting
in basis sets of the :index:`ANO` (:index:`Atomic Natural Orbitals`)
type :cite:`Almlof:87`.

.. compound::

  Some elementary theory: We can do a spectral resolution of a density matrix :math:`D`

  .. math:: D=\sum_k \eta_k c_k c_k^{\text{T}}
     :label: eqn:d-spectral

  where :math:`\eta_k` is the :math:`k`\th eigenvalue (occupation value)
  and :math:`c_k` is the :math:`k`\th eigenvector (natural orbital).
  The occupation number for a natural orbital is a
  measure of how much this orbital contributes to
  the total one-electron density.
  A natural choice is to disregard the natural orbitals
  with small occupation numbers and use those with large
  occupation numbers to form contracted basis functions as

  .. math:: \varphi_k=\sum_i c_{ki} \chi_i

  where :math:`\chi_i` is the :math:`i`\th primitive basis function.

As a generalization to this approach we can
average over density
matrices from several wave functions, resulting
in basis sets of the density matrix averaged ANO type,
see for example :cite:`anoI,anoII,anoIII,anoIV`.
We can view the averaging of density matrices as a sequence
of rank-1 updates in the same way as in equation :eq:`eqn:d-spectral`.
We have more update vectors than the rank of the matrix, but this
does not really change anything. The important observation is
that all :math:`\eta`\s are positive and no information is lost
in the averaging.

The general guideline for which wave functions to include is
based on what you want to be able to describe.
All wave functions you want an accurate description of
should be included in the averaging.

As an example, let us consider the oxygen atom.
We want to be able to describe the atom by itself accurately,
thus a wave function for the atom is needed, usually at the CI level.
In molecular systems, oxygen usually has a negative charge, thus
including :math:`\ce{O-}` is almost mandatory.
A basis set derived from these two wave function is well
balanced for the majority of systems containing oxygen.
A logical conclusion would be that you need to include a few
*molecular* wave functions of systems containing oxygen, but in
practice this is not necessary. This is due to the fact that
the degrees of freedom describing the orbital shape distortion
when forming bonds are virtually identical to the lowest
correlating orbitals.
On the other hand, a few molecular species have oxygen with
positive charge, thus it may be appropriate to include
:math:`\ce{O+}` in the basis set.

.. compound::

  A wide range of specialized basis sets can also be generated,
  for example a molecular basis set describing Rydberg orbitals,
  see the example in the "Tutorials and Examples" part,
  :numref:`TUT:sec:make_rydberg_basis_sets`.
  There is a possibility to create Rydberg orbitals
  automatically by using the keyword
  :kword:`RYDBERG`. Here all unoccupied orbitals with
  negative orbital energies will be used with the associated
  occupation numbers

  .. math:: \eta_k = e^{6.9(\epsilon_k/\epsilon_0-1)}

  where :math:`\epsilon_k` is the orbital energy of orbital :math:`k` and
  :math:`\epsilon_0` is the lowest orbital energy of all
  virtual orbitals. In order to use this option you need
  to use the
  :program:`SCF` or :program:`RASSCF` program to compute
  the orbitals for a cationic system.

You need one or more wave functions,
represented by formatted orbital files,
to generate the average density matrix.
These natural orbital files can be produced by any of the
wave function generators
:program:`SCF`,
:program:`RASSCF`,
:program:`MRCI` or
:program:`CPF`.
You could also use
:program:`MBPT2` or
:program:`CASPT2`.
This approach has been used in the generation of the ANO-RCC basis sets.
Your specific requirements dictate the choice of
wave function generator, but :program:`MRCI` would
be most commonly used.

You are not restricted to atomic calculations but
can mix molecular and atomic calculations freely.
The restrictions are that the name of the center, for which
you are constructing a basis set, must be the same
in all wave functions.
The center may not be "degenerate", i.e.
it may not generate other centers through symmetry
operations. See the description of :program:`SEWARD`
on :numref:`UG:sec:seward`
for a more extensive discussion.
For example for :math:`\ce{O2}` you cannot use :math:`D_{2h}` symmetry
since this would involve one center that is mirrored into the other.
Another restriction is, of course, that you must use the
same primitive set in all calculations.

.. _UG\:sec\:genano_dependencies:

Dependencies
------------

:program:`GENANO` needs one or more wave functions in the
form of natural orbitals. Thus you need to run one or
more of
:program:`SCF`,
:program:`RASSCF`,
:program:`MRCI` or
:program:`CPF`.
You could also use, for example, :program:`MBPT2` or :program:`CASPT2`
but this is in general not recommended.
:program:`GENANO` also needs the one electron file
:file:`ONEINT` and the :file:`RUNFILE` generated by :program:`SEWARD`.

.. index::
   pair: Files; GENANO

.. _UG\:sec\:genano_files:

Files
-----

Below is a list of the files that :program:`GENANO`
reads/writes.
Files :file:`ONEnnn`, :file:`RUNnnn` and :file:`NATnnn` must be supplied to
the program.
Files :file:`ANO` and :file:`FIG` are generated.
File :file:`PROJ` is an optional input file.

Input files
...........

.. class:: filelist

:file:`RUNnnn`
  This file contains miscellaneous information for the nnn'th
  wave function,
  generated by the program :program:`SEWARD`.
  One file per wave function must be supplied,
  :file:`RUN001`, :file:`RUN002`, ....

:file:`ONEnnn`
  This file contains the one-electron integrals corresponding to
  the nnn'th wave function, generated by the program :program:`SEWARD`.
  One file per wave function must be supplied,
  :file:`ONE001`, :file:`ONE002`, ....

:file:`NATnnn`
  This file contains the natural orbitals corresponding to the
  nnn'th wave function, generated by the appropriate wave function
  generating program.
  One file per wave function must be supplied,
  :file:`NAT001`, :file:`NAT002`, ....

:file:`PROJ`
  This file contains orbitals used for projection of the densities.
  Needs to be available if the keyword :kword:`PROJECT`
  is specified.
  It is compatible in format with the file :file:`ANO`, and can thus be the
  the file :file:`ANO` from a previous run of :program:`GENANO`.

Output files
............

.. class:: filelist

:file:`FIG`
  This file contains a PostScript figure file of eigenvalues.

:file:`ANO`
  This file contains the contraction coefficient matrix organized
  such that each column correspond to one contracted basis function.

.. _UG\:sec\:genano_input:

Input
-----

.. compound::

  The input file must contain the line ::

  &GENANO

  right before the actual input starts. Below is a list of the available keywords.
  Please note that you can not abbreviate any keyword.

.. class:: keywordlist

:kword:`TITLE`
  This keyword starts the reading of title lines,
  with no limit on the number of title lines.
  Reading the input as title lines is stopped as soon
  an the input parser detects one of the other keywords.
  This keyword is *optional*.

  .. xmldoc:: <KEYWORD MODULE="GENANO" NAME="TITLE" APPEAR="Title" LEVEL="BASIC" KIND="STRING">
              %%Keyword: TITLe <basic>
              <HELP>
              This keyword starts the reading of title lines,
              with no limit on the number of title lines.
              Reading the input as title lines is stopped as soon
              an the input parser detects one of the other keywords.
              </HELP>
              This keyword is optional.
              </KEYWORD>

:kword:`SETS`
  This keyword indicates that the next line of input
  contains the number of sets to be used in the
  averaging procedure.
  This keyword must precede :kword:`WEIGHTS` if
  both are supplied.
  This keyword is *optional*, with one set as the default.

  .. xmldoc:: <KEYWORD MODULE="GENANO" NAME="SETS" APPEAR="Sets" LEVEL="BASIC" KIND="INT">
              %%Keyword: SETS <basic>
              <HELP>
              This keyword indicates that the next line of input
              contains the number of sets to be used in the
              averaging procedure.
              </HELP>
              This keyword must precede keyword WEIGHTS if
              both are supplied.
              This keyword is optional, with one set as the default.
              </KEYWORD>

:kword:`CENTER`
  This keyword is followed, on the next line, by the atom
  label for which the basis set is to be generated.
  The label must match the label you supplied to
  :program:`SEWARD`.
  In previous versions of :program:`GENANO` this label had to
  be in uppercase, but this restriction is now lifted and
  the case does not matter.
  This keyword is *compulsory*.

  .. xmldoc:: <KEYWORD MODULE="GENANO" NAME="CENTER" APPEAR="Center" LEVEL="BASIC" KIND="STRING">
              %%Keyword: CENTer <basic>
              <HELP>
              This keyword is followed, on the next line, by the atom
              label for which the basis set is to be generated.
              The label must match the label you supplied to
              SEWARD.
              </HELP>
              In previous versions of GENANO this label had to
              be in uppercase, but this restriction is now lifted and
              the case does not matter.
              This keyword is compulsory.
              </KEYWORD>

:kword:`ROWWISE`
  This keyword makes :program:`GENANO` produce the
  contraction coefficients row-wise instead of
  column-wise as is the default.
  This keyword is *optional*.

  .. xmldoc:: <KEYWORD MODULE="GENANO" NAME="ROWWISE" APPEAR="Row-wise" LEVEL="BASIC" KIND="SINGLE">
              %%Keyword: ROWWise <advanced>
              <HELP>
              This keyword makes GENANO to produce the
              contraction coefficients row-wise instead of
              column-wise as is the default.
              </HELP>
              This keyword is optional.
              </KEYWORD>

:kword:`WEIGHTS`
  This keyword must be subsequent to keyword :kword:`SETS`
  if both are supplied.
  This keyword is *optional*,
  with equal weight on each of the sets as default.

  .. xmldoc:: <KEYWORD MODULE="GENANO" NAME="WEIGHTS" APPEAR="Weights" LEVEL="BASIC" KIND="REALS_LOOKUP" SIZE="SETS">
              %%Keyword: WEIGhts <basic>
              <HELP>
              </HELP>
              This keyword must be subsequent to keyword SETS
              if both are supplied.
              This keyword is optional,
              with equal weight on each of the sets as default.
              </KEYWORD>

:kword:`PROJECT`
  This keyword states that you want to project out certain
  degrees of freedom from the density matrix.
  This can be useful for generating, for example,
  node less valence orbitals to be used with ECP's.
  If this keyword is specified, you must supply the file
  :file:`PROJ` obtained as file :file:`ANO` from a previous
  :program:`GENANO` calculation, for instance.
  This keyword is *optional*.

  .. xmldoc:: <KEYWORD MODULE="GENANO" NAME="PROJECT" APPEAR="Project out" LEVEL="ADVANCED" KIND="SINGLE">
              %%Keyword: PROJect <advanced>
              <HELP>
              This keyword states that you want to project out certain
              degrees of freedom from the density matrix.
              This can be useful for generating, for example,
              nodeless valence orbitals to be used with ECP's.
              If this keyword is specified, you must supply the file
              PROJ obtained as file ANO from a previous
              GENANO calculation, for instance.
              </HELP>
              This keyword is optional.
              </KEYWORD>

:kword:`LIFTDEGENERACY`
  This keyword will modify the occupation numbers read from
  the orbitals files. The purpose is to lift the
  degeneracy of core orbitals to avoid rotations.
  The occupation numbers are changed according to
  :math:`\eta'=\eta(1+10^{-3}/n)`
  where :math:`n` is the sequence number of the orbital
  in its irreducible representation.
  This keyword is *optional*.

  .. xmldoc:: <KEYWORD MODULE="GENANO" NAME="LIFTDEGENERACY" APPEAR="Lift degeneracy" LEVEL="ADVANCED" KIND="SINGLE">
              %%Keyword: LIFTdegeneracy <advanced>
              <HELP>
              This keyword will modify the occupation numbers read from
              the orbitals files. The purpose is to lift the
              degeneracy of core orbitals to avoid rotations.
              The occupation numbers are changed according to
              o'=o*(1+10^-3/n)
              where n is the sequence number of the orbital
              in its irreducible representation.
              </HELP>
              This keyword is optional.
              </KEYWORD>

:kword:`RYDBERG`
  This keyword enables automatic generation of Rydberg
  orbitals. With this keyword all occupied orbitals
  will get occupation number zero while the virtual
  orbitals will get a small occupation number
  decreasing with orbital number. Useful with a calculation
  on an cation where the virtual orbitals are near perfect
  Rydberg orbitals.
  Note that you must use orbitals from the
  :program:`SCF` or
  :program:`RASSCF` program.
  This keyword is *optional*.

  .. xmldoc:: <KEYWORD MODULE="GENANO" NAME="RYDBERG" APPEAR="Rydberg orbitals" LEVEL="ADVANCED" KIND="SINGLE">
              %%Keyword: RYDBerg <advanced>
              <HELP>
              This keyword enables automatic generation of Rydberg orbitals.
              With this keyword all occupied orbitals will get occupation
              number zero while the virtual orbitals will get a small
              occupation number decreasing with orbital number. Useful
              with a calculation on an cation where the virtual orbitals
              are near perfect Rydberg orbitals. Note that you must use
              orbitals from the SCF or RASSCF program.
              </HELP>
              This keyword is optional.
              </KEYWORD>

:kword:`NOTHRESHOLD`
  This keyword is used to specify the threshold for
  keeping NO's (natural orbitals). Orbitals with
  occupation numbers less than the threshold are
  discarded. The threshold is read from the line
  following the keyword. Default value is 1.0d-8.

  .. xmldoc:: <KEYWORD MODULE="GENANO" NAME="NOTHRESHOLD" APPEAR="Natural orbital threshold" LEVEL="ADVANCED" KIND="REAL" DEFAULT_VALUE="1.0d-8">
              %%Keyword: NOTHreshold <advanced>
              <HELP>
              This keyword is used to specify the threshold for
              keeping NO's (natural orbitals). Orbitals with
              occupation numbers less than the threshold are
              discarded. The threshold is read from the line
              following the keyword.
              </HELP>
              Default value is 1.0d-8.
              </KEYWORD>

Below is a simple input example, where we construct an
ANO basis set for the carbon atom.
Two wave functions are used, the SCF wave function and the
SDCI wave function for the ground state of the atom.

.. extractfile:: ug/GENANO.input

  &SEWARD
  Title
   Carbon atom
  Symmetry
  x y z
  Expert
  Basis set
  C..... / inline
    6.0 2
     10   10
  5240.6353 782.20479 178.35083 50.815942 16.823562 6.1757760 2.4180490
  .51190000 .15659000 .05480600
  1. 0. 0. 0. 0. 0. 0. 0. 0. 0.
  0. 1. 0. 0. 0. 0. 0. 0. 0. 0.
  0. 0. 1. 0. 0. 0. 0. 0. 0. 0.
  0. 0. 0. 1. 0. 0. 0. 0. 0. 0.
  0. 0. 0. 0. 1. 0. 0. 0. 0. 0.
  0. 0. 0. 0. 0. 1. 0. 0. 0. 0.
  0. 0. 0. 0. 0. 0. 1. 0. 0. 0.
  0. 0. 0. 0. 0. 0. 0. 1. 0. 0.
  0. 0. 0. 0. 0. 0. 0. 0. 1. 0.
  0. 0. 0. 0. 0. 0. 0. 0. 0. 1.
      6    6
  18.841800 4.1592400 1.2067100 .38554000 .12194000 .04267900
  1. 0. 0. 0. 0. 0.
  0. 1. 0. 0. 0. 0.
  0. 0. 1. 0. 0. 0.
  0. 0. 0. 1. 0. 0.
  0. 0. 0. 0. 1. 0.
  0. 0. 0. 0. 0. 1.
      3    3
  1.2838000 .34400000 .09220000
  1. 0. 0.
  0. 1. 0.
  0. 0. 1.
  C  0.000000  0.000000  0.000000
  End of basis

  &SCF
  Occupied =  2 0 0 0 0 0 0 0

  &RASSCF
  Symmetry =  4
  Spin     =  3
  nActEl   =  2 0 0
  Frozen   =  0 0 0 0 0 0 0 0
  Inactive =  2 0 0 0 0 0 0 0
  Ras2     =  0 1 1 0 0 0 0 0
  LevShft  =  0.00
  LumOrb
  Thrs     =  0.1d-8 0.1d-4 0.1d-4

  &MOTRA
  LumOrb
  Frozen   =  1 0 0 0 0 0 0 0

  &GUGA
  Electrons =  4
  Spin      =  3
  Inactive  =  1 0 0 0 0 0 0 0
  Active    =  0 1 1 0 0 0 0 0
  CiAll     =  4

  &MRCI
  SDCI

  >>COPY $Project.RunFile RUN001
  >>COPY $Project.RunFile RUN002
  >>COPY $Project.OneInt  ONE001
  >>COPY $Project.OneInt  ONE002
  >>COPY $Project.RasOrb  NAT001
  >>COPY $Project.CiOrb   NAT002

  &GENANO
  Title
   Carbon atom
  Project
  sets
   2
  Center
  C
  Weights
   0.5 0.5
  >>RM ONE001
  >>RM ONE002
  >>RM NAT001
  >>RM NAT002

.. xmldoc:: </MODULE>
