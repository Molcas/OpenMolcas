.. index::
   single: Program; RASSI
   single: RASSI

.. _UG\:sec\:rassi:

:program:`rassi`
================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="RASSI">
            %%Description:
            <HELP>
            %%Description:
            The RASSI program calculates overlaps, and matrix
            elements of one-electron operators, and of the electronic Hamiltonian,
            over a basis of RASSCF wave functions, which may each have its own
            independent set of orbitals. Energies and matrix elements are
            computed also for the non-interacting linear combinations of states,
            i.e., doing a limited CI using the RASSCF states as a non-orthogonal basis.
            RASSI is extensively used for computing dipole oscillator strengths.
            Finally, it can also compute e.g. spin-orbit interaction matrix elements,
            transition dipole moments, (bi-)natural orbitals and other quantities.
            </HELP>

The
:program:`RASSI` (RAS :index:`State Interaction`) program forms overlaps and
other :index:`matrix
elements <single: Matrix elements; RASSCF>` of the Hamiltonian and other operators
over a wave function basis, which consists of RASSCF wave functions,
each with an individual set of orbitals. It is extensively used
for computing dipole :index:`oscillator strengths <single: Oscillator strength>`, but any
one-electron operator, for which the :program:`Seward` has computed
integrals to the :file:`ORDINT` file, can be used, not just dipole
moment components.

Also, it solves the Schrödinger
equation projected on the space spanned by these wave functions,
i.e., it forms non-interacting linear combinations of the input
state functions, and computes matrix elements over the resulting
eigenbasis as well.

Finally, using these spin-free eigenstates as a basis, it can
compute spin--orbit interaction matrix elements, diagonalize
the resulting matrix, and compute various matrix elements over
the resulting set of spin--orbit eigenstates.

If only matrix
elements of some one-electron operator(s), such as the dipole
transition moments, are required, the calculation of Hamiltonian
matrix elements and the transformation to the eigenbasis of this
matrix can be skipped. However, if any states have the same symmetry
and different orbitals, it is desirable to use the transitions strengths
as computed between properly non-interacting and orthonormal states.
The reason is that the individually optimized RASSCF states are
interacting and non-orthogonal, and the main error in the computed
transition matrix elements is the difference in electronic dipole
moment times the overlap of any two states involved. For excited
states, the overlap is often in the order of 10%.

Please note: Due to the increasing number of calculations done with
a hundred input states, or more, there has been a demand to change
the output. Until |molcas| 6.2, the default assumption has been to print
all expectation values and matrix elements that can be computed from
the selection of one-electron integrals. From 6.4, this is requested by
keywords, see the keyword list below for XVIN, XVES, XVSO, MEIN,
MEES, and MESO.

Apart from computing oscillator strengths, overlaps and Hamiltonian
matrix elements can be used to compute :index:`electron transfer rates <single: Electron transfer rate>`, or
to form :index:`quasi-diabatic states <single: Quasi-diabatic states>` and reexpress matrix elements over a
basis of such states.

The CSF space of a RASSCF wave function is closed under deexcitation.
For any given pair of RASSCF wave functions, this is used in the
way described in reference :cite:`Malmqvist:86` to allow the pair of orbital
sets to be transformed to a biorthonormal pair, while simultaneously
transforming the CI expansion coefficients so that the wave functions
remain unchanged. The basic principles are the same as in the earlier
program :cite:`Malmqvist:89`, but is adapted to allow RASSCF as well as
CASSCF wave functions. It uses internally a Slater determinant
expansion. It can now use spin-dependent operators,
including the AMFI spin--orbit operator, and can compute matrix elements
over spin--orbit states, i.e. the eigenstates of the sum of the
spin-free hamiltonian and the spin--orbit operator.

One use of the RASSI eigenstates is to resolve ambiguities due
to the imperfect description of highly excited states.
Association between individually optimized states and the exact
electronic eigenstates is often not clear, when the calculation
involves several or many excited states. The reason is that the
different states each use a different set of orbitals. The State
Interaction calculation gives an unambiguous set of non-interacting and
orthonormal eigenstates to the projected Schrödinger equation, and
also the overlaps between the original RASSCF wave functions and the
eigenstates. The latter is a very efficient diagnostic, since it
describes the RASSCF states in terms of one single wave-function basis
set.

.. index::
   single: RASSCF; Multiple solutions

To make the last point clear, assume the following situation:
We have
performed three RASSCF calculations, one where we optimize for the
lowest state, one for the first excited state, and one for the 2nd
excited state in the same symmetry. The active orbitals are fairly
much mixed around, so a simple inspection of the CI coefficient is
insufficient for comparing the states. Assume that for each state, we
have calculated the three lowest CI roots. It can now happen, that the
2nd root of each calculation is a fair approximation to the exact 2nd
eigenstate, and the same with the 3rd, or possibly that the order gets
interchanged in one or two of the calculation. In that case, a RASSI
calculation with these 9 states will give three improved solutions
close to the original ones, and of course 6 more that are considered
to be the removed garbage. The overlaps will confirm that each of the
input states consists mainly of one particular out of the three lowest
eigenstates. This situation is the one we usually assume, if no
further information is available.

However, it happens that the active orbitals of the three calculations
do not span approximately the same space. The orbital optimization
procedure has made a qualitatively different selection of correlating
orbitals for the three different calculation. Then the :program:`RASSI`
calculation may well come out with 4 lowest roots that overlap
strongly with the original RASSCF states. This may change the
assignments and may also give valuable information about the
importance of some state. The natural orbitals of the eigenstates will
show that the active space used in the RASSCF was in some way
inappropriate.

Another bothersome situation is also solved by the RASSI method. The
analysis of the original states in terms of RASSI eigenstates may
reveal that the three optimized RASSCF states consists mainly of TWO
low RASSI eigenstates! This is because the RASSCF optimization
equations are non-linear and may sometimes offer spurious extra
solutions. Two of the calculations are in this case to be regarded
qualitatively, as two different (local) solutions that
approximate (imperfectly) the same excited state. Also in this case, the
natural orbitals will probably offer a clue to how to get rid of the
problem. Extra solutions rarely occur for low states in CASSCF
calculations, provided a generous active space can be afforded.
Problems occur when the active space is too small, and in
particular with general RASSCF calculations.

A further application is the preparation of a suitable orbital basis
for a subsequent CI calculation. Note that such an application also
allows the use of badly converged RASSCF wave functions, or of RASSCF
wave functions containing multiple minima solutions close to a common
exact eigenstate. In effect, the :program:`RASSI` program cleans up the situation
by removing the errors due to bad convergence (pushing the errors into
a garbage part of the spectrum). This requires that the set of input
states (9 in this example) provides flexibility enough to remove at
least a major part of the error. As one would expect, this is usually
true: The erratic non-convergent, or the too slowly convergent, error
mode is to a large extent spanned by the few lowest RASSCF wave
functions.

Finally, there are situations where there is no problem to obtain
adiabatic RASSCF solutions, but where it is still imperative to use
RASSI natural orbitals in a subsequent CI. Consider the case of
transition metal chemistry, where there is in general two or more
electronic states involved. These states are supposed to interact
strongly, at least within some range of interatomic distances. Here,
an MCSCF solution, such as RASSCF, will have at least two very
different solutions, one associated with each configuration of the
transition metal atom. Using one set of orbitals, one electronic state
has a reasonably described potential energy curve, while other states
get pushed far up in energy. Using another set of orbitals, another
state gets correctly described. In no calculation with a single
orbital set do we obtain the avoided crossings, where one switches
from one diabatic state to another. The only way to accomplish this is
via a RASSI calculation. In this case, it is probably necessary also to
shift the energies of the RASSCF states to ensure that the crossing
occur at the correct places. The shifts can be determined by
correcting the atomic spectrum in the separated-atoms limit.

Note, however, that most of the problems described above can be
solved by performing state-averaged RASSCF calculations.

.. In the latest version of |molcas|, derivatives of transition dipole moments
   have been added to RASSI :cite:`Bernhardsson:99b`. The derivatives are
   calculated as the matrix element of the product of the (CI/orbital)
   rotation operator and the dipole moment operator.

.. index::
   pair: Dependencies; RASSI

.. _UG\:sec\:rassi_dependencies:

Dependencies
------------

The :program:`RASSI` program needs one or more :file:`JOBIPH` files produced
by the :program:`RASSCF` program. Also, it needs a :file:`ONEINT` file from
:program:`SEWARD`, with overlap integrals and any one-electron
property integrals for the requested matrix elements. If Hamiltonian
matrix elements are used, also the :file:`ORDINT` file is needed.

.. For derivatives the :program:`RASSI` needs the :file:`MCKINT` file
   produced by :program:`MCKINLEY` and :program:`MCLR` containing
   the derivatives of the dipole operator and the orbital rotations and the
   state transfer operators.

   The existence of a file with the name :file:`MCKINT1` will automatically
   change :program:`RASSI` to from ordinary calculation of state interaction
   properties, to calculation of derivatives of state interaction properties,
   like transition dipole derivatives and non adiabatic coupling constants.

   It is important that if derivatives of state interaction properties should
   be calculated, the expansion center for that property must be (0,0,0).
   For derivatives of transition dipole moments, the following keyword has to be
   added to the :program:`SEWARD` input. ::

     Center= 1; 1 0.0 0.0 0.0

.. index::
   pair: Files; RASSI

.. _UG\:sec\:rassi_files:

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

:file:`JOBnnn`
  A number of :file:`JOBIPH` files from different :program:`RASSCF` jobs.
  An older naming convention assumes file names JOB001, JOB002, etc. for these files.
  They are automatically linked to default files named :file:`$Project.JobIph`,
  :file:`$Project.JobIph01`, :file:`$Project.JobIph02`, etc. in directory :file:`$WorkDir`,
  unless they already exist as files or links before the program starts.
  You can set up such links yourself, or else you can specify file names
  to use by the keyword :kword:`IPHNames`.

:file:`JOBIPHnn`
  A number of :file:`JOBIPH` files from different :program:`RASSCF` jobs.
  The present naming convention assumes file names JOBIPH, JOBIPH01, etc. for
  such files, when created by subsequent :program:`RASSCF` runs, unless
  other names were specified by input.
  They are automatically linked to default files named :file:`$Project.JobIph`,
  :file:`$Project.JobIph01`, :file:`$Project.JobIph02`, etc. in directory :file:`$WorkDir`,
  unless they already exist as files or links before the program starts.
  You can set up such links yourself, or else you can specify file names
  to use by the keyword :kword:`IPHNames`.

  .. :file:`MCKINTn`
       A number of :file:`MCKINT` files from different Single state calculations.
       The numbering of the MCKINTn files should be the same as the
       numbering of JOBnnn files, files with the same number should correspond to
       the same state.

Output files
............

.. class:: filelist

:file:`SIORBnn`
  A number of files containing natural orbitals, (numbered sequentially as
  :file:`SIORB01`, :file:`SIORB02`, etc.)

:file:`BIORBnnmm`
  A number of files containing binatural orbitals for the transition between
  states ``nn`` and ``mm``. Each such file contains pairs of orbitals, in the same format
  as the :math:`\alpha` and :math:`\beta` components of UHF orbitals. The file for transition
  to state ``nn``\ =2 from state ``mm``\ =1 will be named :file:`BIORB.2_1`.

:file:`TOFILE`
  This output is only created if :kword:`TOFIle` is given in the input.
  It will contain the transition density matrix computed by :program:`Rassi`.
  Currently, this file is only used as input to :program:`QmStat`.

:file:`EIGV`
  Like :file:`TOFILE` this file is only created if :kword:`TOFIle` is given
  in the input. It contains auxiliary information that is picked up
  by :program:`QmStat`.

  .. :file:`UNSYM`
       The derivative of the transition dipole moment desymmetrized.

.. index::
   pair: Input; RASSI

.. _UG\:sec\:rassi_input:

Input
-----

This section describes the input to the
:program:`RASSI` program in the |molcas| program system,
with the program name: ::

  &RASSI

When a keyword is followed by additional mandatory lines of input,
this sequence cannot be interrupted by a comment line. The first 4
characters of keywords are decoded. An unidentified keyword makes the
program stop.

.. index::
   pair: Keywords; RASSI

Keywords
........

.. :kword:`CHOLesky`
     :program:`RASSI` will use Cholesky (or RI/DF) representation of the two-electron integrals to compute
     the corresponding contributions to the Fock matrices and to perform the MO integral transformations.
     The default (LK) algorithm is used. The configuration may be tailored using the ChoInput section.
     Default is to not use Cholesky unless the Cholesky (or RI/DF) representation of the two-electron
     integrals has been produced by :program:`SEWARD`.

..   .. xmldoc:: %%Keyword: Cholesky <basic>
                 Use of Cholesky (or RI/DF) representation for the two-electron integrals
                 with default RASSI settings.

.. class:: keywordlist

:kword:`CHOInput`
  This marks the start of an input section for modifying
  the default settings of the Cholesky RASSI.
  Below follows a description of the associated options.
  The options may be given in any order,
  and they are all optional except for
  :kword:`ENDChoinput` which marks the end of the :kword:`CHOInput` section.

  .. xmldoc:: <GROUP MODULE="RASSI" NAME="CHOINPUT" APPEAR="Cholesky input section" KIND="BLOCK" LEVEL="ADVANCED">
              %%Keyword: Choinput <advanced>
              <HELP>
              Manually modify the settings of the Cholesky RASSI.
              </HELP>

  * :kword:`NoLK`
    Available only within ChoInput. Deactivates the "Local Exchange" (LK) screening algorithm :cite:`Aquilante:07a` in computing
    the Fock matrix. The loss of speed compared to the default algorithm can be substantial, especially for electron-rich systems.
    Default is to use LK.

    .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="NOLK" APPEAR="Turn Off LK screening" LEVEL="ADVANCED" KIND="SINGLE" REQUIRE="CHOI">
                %%Keyword: NoLK <advanced>
                <HELP>
                Deactivates LK screening.
                </HELP>
                </KEYWORD>

  * :kword:`DMPK`
    Available only within ChoInput. Modifies the thresholds used in the LK screening.
    The keyword takes as argument a (double precision) floating point (non-negative) number used
    as correction factor for the LK screening thresholds.
    The default value is 1.0d-1. A smaller value results in a slower but more accurate calculation.

    **Note:** the default choice of the LK screening thresholds is tailored to achieve as much as possible an
    accuracy of the RASSI energies consistent with the choice of the Cholesky decomposition
    threshold.

    .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="DMPK" APPEAR="Damping for LK" LEVEL="ADVANCED" KIND="REAL" EXCLUSIVE="NOLK" REQUIRE="CHOI">
                %%Keyword: DMPK <advanced>
                <HELP>
                Modifies the thresholds used in the LK screening.
                The default value is 1.0d-1. A smaller value results in a slower but more accurate calculation.
                </HELP>
                </KEYWORD>

  * :kword:`NODEcomposition`
    Available only within ChoInput. The inactive Exchange contribution to the Fock matrix is computed using inactive canonical orbitals
    instead of (localized) "Cholesky MOs".
    This choice is effective only in combination with the LK screening.
    Default is to use Cholesky MOs. **Note:** the Cholesky MOs in RASSI are computed by decomposing the
    density type supermatrix :math:`\mat{D}=(\mat{C}_A, \mat{C}_B)(\mat{C}_A, \mat{C}_B)^{\text{T}}` where :math:`\mat{C}` is the corresponding canonical
    MOs matrix for the state :math:`A` and :math:`B`.

    .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="NODE" APPEAR="Turn Off density decomposition" LEVEL="ADVANCED" KIND="SINGLE" REQUIRE="CHOI">
                %%Keyword: NODE <advanced>
                <HELP>
                The inactive Exchange contribution to the Fock matrix is computed using inactive canonical orbitals
                instead of (localized) "Cholesky MOs".
                </HELP>
                </KEYWORD>

  * :kword:`PSEUdo`
    When computing the coupling between 2 different
    states A and B, only for the first state we use pure Cholesky MOs. The invariance of the Fock matrix
    is then ensured by rotating the orbitals of B according to the orthogonal matrix defined in A
    through the Cholesky localization. These orbitals used for B are therefore called "pseudo Cholesky MOs".

    .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="PSEU" APPEAR="Use PseudoCholesky orbitals" LEVEL="ADVANCED" KIND="SINGLE" REQUIRE="CHOI">
                %%Keyword: PSEU <advanced>
                <HELP>
                The inactive Exchange contribution to the Fock matrix is computed using pseudo Cholesky orbitals.
                </HELP>
                </KEYWORD>

    .. xmldoc:: </GROUP>

  * :kword:`TIME`
    Activates printing of the timings of each task of the Fock matrix build.
    Default is to not show these timings.

  * :kword:`MEMFraction`
    Set the fraction of memory to use as global Cholesky vector buffer.
    Default: for serial runs 0.0d0; for parallel runs 0.3d0.

:kword:`MEIN`
  Demand for printing matrix elements of all selected one-electron
  properties, over the input RASSCF wave functions.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="MEIN" APPEAR="RASSCF matrix elements" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: MEIN <basic>
              <HELP>
              Demand for printing matrix elements of all selected one-electron
              properties, over the input RASSCF wave functions.
              </HELP>
              </KEYWORD>

:kword:`MEES`
  Demand for printing matrix elements of all selected one-electron
  properties, over the spin-free eigenstates.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="MEES" APPEAR="spin-free matrix elements" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: MEES <basic>
              <HELP>
              Demand for printing matrix elements of all selected one-electron
              properties, over the spin-free eigenstates.
              </HELP>
              </KEYWORD>

:kword:`MESO`
  Demand for printing matrix elements of all selected one-electron
  properties, over the spin--orbit states.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="MESO" APPEAR="SO matrix elements" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: MESO <basic>
              <HELP>
              Demand for printing matrix elements of all selected one-electron
              properties, over the spin-orbit states.
              </HELP>
              </KEYWORD>

  .. :kword:`PRINT`
     Set individual print levels for various subroutines of the code, mainly
     for debugging purposes.
     This keyword requires an entry with number of name,value
     pairs, with the same format as that used for the PROP input. The names
     are subroutine names, and each value is the print level setting for that
     subroutine.

  .. .. xmldoc:: %%Keyword: Print <advanced>
                 Set individual print levels for various subroutines of the code.

:kword:`PROPerty`
  Replace the default selection of one-electron operators, for which
  matrix elements and expectation values are to be calculated, with a
  user-supplied list of operators.

  .. compound::

    From the lines following the keyword the selection list is
    read by the following *FORTRAN* code: ::

      READ({*},{*}) NPROP,(PNAME(I),ICOMP(I),I=1,NPROP)

    NPROP is the number of selected properties, PNAME(I) is a
    character string with the label of this operator on :program:`SEWARD`'s
    one-electron integral file, and ICOMP(I) is the component number.

  The default selection is to use dipole and/or velocity integrals, if
  these are available in the :file:`ONEINT` file. This choice is replaced by the
  user-specified choice if the :kword:`PROP` keyword is used.
  Note that the character strings are read using list directed input and
  thus must be within single quotes, see sample input below.
  For a listing of presently available operators, their labels, and
  component conventions, see
  :program:`SEWARD` program description.

  .. xmldoc:: %%Keyword: Property <basic>
              Enter a user-supplied selection of one-electron operators, for which
              matrix elements and expectation values are to be calculated. Without
              this keyword and list, the default choice is to use every operator,
              for which the one-electron integral file supplies integrals.

:kword:`SOCOupling`
  Enter a positive threshold value. Spin--orbit interaction matrix
  elements over the spin components of the spin-free eigenstates
  will be printed, unless smaller than this threshold.
  The value is given in cm\ :math:`^{-1}` units. The keyword is
  ignored unless an SO hamiltonian is actually computed.

  .. xmldoc:: %%Keyword: SOCoupling <basic>
              Enter a positive threshold value. Spin-orbit interaction matrix
              elements over the spin components of the spin-free eigenstates
              will be printed, unless smaller than this threshold.
              The value is given in cm-1 units. The keyword is ignored unless
              an SO hamiltonian is actually computed.

:kword:`SOPRoperty`
  Enter a user-supplied selection of one-electron operators, for which
  matrix elements and expectation values are to be calculated over the
  spin--orbit eigenstates. This keyword has no effect unless the
  :kword:`SPIN` keyword has been used. Format: see :kword:`PROP` keyword.

  .. xmldoc:: %%Keyword: SOProperty <basic>
              Enter a selection of one-electron operators, for which
              matrix elements and expectation values are to be calculated over the
              spin-orbit eigenstates. This keyword has no effect unless the
              SPIN keyword has been used. Format: see PROP keyword.

:kword:`SPINorbit`
  Spin--orbit interaction matrix elements will be computed. Provided that
  the :kword:`ONEL` keyword was not used, the resulting Hamiltonian including the
  spin--orbit coupling, over a basis consisting of all the spin components
  of wave functions constructed using the spin-free eigenstates, will be
  diagonalized. NB: For this keyword to have any effect, the SO integrals
  must have been computed by :program:`SEWARD`! See :kword:`AMFI` keyword in :program:`SEWARD` documentation.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="SPIN" APPEAR="spin-orbit calc." KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: Spinorbit <basic>
              <HELP>
              Spin-orbit interaction matrix elements will be computed and the resulting
              Hamiltonian including will be diagonalized.
              NB: For this keyword to have any effect, the SO integrals
              must have been computed by SEWARD (AMFI integrals)!
              </HELP>
              </KEYWORD>

:kword:`ONEL`
  The two-electron integral file will not be accessed. No Hamiltonian
  matrix elements will be calculated, and only matrix elements for the
  original RASSCF wave functions will be calculated.
  :kword:`ONEE` is a valid synonym for this keyword.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="ONEL" APPEAR="One-electron properties only" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: Onel <basic>
              <HELP>
              The two-electron integral file will not be accessed. No hamiltonian
              matrix elements will be calculated. Non-interacting states will not
              be formed. Onee is a valid synonym for this keyword.
              </HELP>
              </KEYWORD>

:kword:`J-VAlue`
  For spin--orbit calculations with single atoms, only: The output lines
  with energy for each spin--orbit state will be annotated with the
  approximate J and Omega quantum numbers.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="J-VALUE" APPEAR="J-Value" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: J-Value <basic>
              <HELP>
              For spin-orbit calculations with single atoms, only: The output lines
              with energy for each spin-orbit state will be annotated with the
              approximate J quantum number.
              </HELP>
              </KEYWORD>

:kword:`OMEGa`
  For spin--orbit calculations with linear molecules, only: The output lines
  with energy for each spin--orbit state will be annotated with the
  approximate Omega quantum number.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="OMEGA" APPEAR="Omega-Value" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: Omega <basic>
              <HELP>
              For spin-orbit calculations on linear molecules, only: The output lines
              with energy for each spin-orbit state will be annotated with the
              approximate Omega quantum number.
              </HELP>
              </KEYWORD>

:kword:`NROF jobiphs`
  Number of
  :file:`JOBIPH` files used as input. This keyword should be
  followed by the number of
  states to be read from each :file:`JOBIPH`. Further, one line per
  :file:`JOBIPH` is required with a list of the states to be
  read from the particular file. See sample input below.
  Alternatively, the first line can contain the number of :file:`JOBIPH` used
  as input followed by the word "``ALL``", indicating that all states
  will be taken from each file. In this case no further lines are required.
  For :file:`JOBIPH` file names, see the Files section.
  Note: If this keyword is missing, then by default all files named "JOB001",
  "JOB002", etc. will be used, and all states found on these files will be
  used.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="JOBIPH" APPEAR="Input states from JOBIPHs" KIND="STRINGS" LEVEL="BASIC">
              %%Keyword: NrOf <basic>
              <HELP>
              Number of JOBIPH files used as input, followed by a list of
              the number of states to read from each JOBIPH, and finally,
              for each JOBIPH a list of which states to select.
              </HELP>
              </KEYWORD>


:kword:`REDL`
  In many cases, RASSI is used to compute the transition moments between
  a set of initial states (for example the ground state) and a set of final states.
  "Reduced loop" allows to restrict the computation of transition moments between the two sets
  and not within each set, thus saving time and reducing the output size.
  The keyword is followed by the index where the two sets split (assuming energy ordering).
  For a calculation between one ground state and several excited states, REDL should be 1.
  Default is to compute the transition moments between all states.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="REDL" APPEAR="Reduced loop" KIND="INT" SIZE="1" LEVEL="BASIC">
              %%Keyword: Redl <basic>
              <HELP>
              Restricts the computation of transition moments to be only between
              two sets of states, and not also within each set.
              The keyword is followed by the number of states
              in the first set (assuming energy ordering).
              </HELP>
              </KEYWORD>

:kword:`IPHNames`
  Followed by one entry for each :file:`JOBIPH` file to be used, with the
  name of each file. Note: This keyword presumes that the number of
  :file:`JOBIPH` files have already been entered using keyword :kword:`NROF`.
  For default :file:`JOBIPH` file names, see the Files section.
  The names will be truncated to 8 characters and converted to uppercase.

  .. xmldoc:: %%Keyword: IPHNames <basic>
              Followed by one entry for each JOBIPH file to be used, with the
              name of each file. Note: This keyword presumes that the number of
              JOBIPH files have already been entered using keyword NROF.
              The names will be truncated to 8 characters and converted to uppercase.

  .. :kword:`NACMe`
       Switch from calculations of transition dipole moment to calculation of
       non adiabatic coupling constants. Only valid if :file:`MCKINT1` file exist.

       .. .. xmldoc:: %%Keyword: NACM <advanced>
                      For computing non-adiabatic coupling matrix elements. This requires
                      that a MCKINT1 file exist. After this keyword should follow a list
                      of energy shifts, one for each wave function. Such shifts are
                      usually needed in order to ensure that energy crossings occur where
                      they should. Note: this keyword must not precede the NROF input.

:kword:`SHIFt`
  The next entry or entries gives an energy shift for each wave function,
  to be added to diagonal elements of the Hamiltonian matrix.
  This may be necessary e.g. to ensure that an energy crossing occurs
  where it should. NOTE: The number of states must be known
  (See keyword :kword:`NROF`) before this input is read.
  In case the states are not orthonormal, the actual quantity added to
  the Hamiltonian is ``0.5D0*(ESHFT(I)+ESHFT(J))*OVLP(I,J)``. This is necessary
  to ensure that the shift does not introduce artificial interactions.
  :kword:`SHIFT` and :kword:`HDIAG` can be used together.

  .. xmldoc:: %%Keyword: Shift <basic>
              The next entry or entries gives an energy shift for each wave function,
              to be added to diagonal elements of the Hamiltonian matrix.

:kword:`HDIAg`
  The next entry or entries gives an energy for each wave function,
  to replace the diagonal elements of the Hamiltonian matrix.
  Non-orthogonality is handled similarly as for the :kword:`SHIFT` keyword.
  :kword:`SHIFT` and :kword:`HDIAG` can be used together.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="HDIAG" APPEAR="Diagonal elements" KIND="REALS_COMPUTED" SIZE="1" LEVEL="BASIC">
              %%Keyword: HDiag <basic>
              <HELP>
              Enter an energy for each spin-free wave function,
              to replace the diagonal elements of the Hamiltonian matrix.
              For example to use CASPT2 shifted energies in the diagonal.
              </HELP>
              </KEYWORD>

:kword:`NATOrb`
  The next entry gives the number of eigenstates for which natural
  orbitals will be computed. They will be written, formatted, commented,
  and followed by natural occupancy numbers, on one file each state.
  For file names, see the Files section.
  The format allows their use as standard orbital input files to
  other |molcas| programs.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="NATORB" APPEAR="Natural Orbitals" KIND="INT" LEVEL="BASIC">
              %%Keyword: NatOrb <basic>
              <HELP>
              Enter the number of eigenstates for which natural orbitals should
              be computed and written to file. These will be written together with
              natural occupation numbers in the usual format used by MOLCAS.
              </HELP>
              </KEYWORD>

:kword:`BINAtorb`
  The next entry gives the number of transitions for which binatural
  orbitals will be computed. Then a line should follow for each transition,
  with the two states involved. The orbitals and singular values provide a
  singular value decomposition of a transition density matrix \cite{Malmqvist:2012}.
  The bra and ket orbitals are written followed by the singular values in the
  usual UHF format used by other |molcas| programs.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="BINATORB" APPEAR="Binatural Orbitals" KIND="INTS_COMPUTED" SIZE="2" LEVEL="BASIC">
              %%Keyword: BiNatOrb <basic>
              <HELP>
              Enter the number of transitions, for which binatural orbitals should
              be computed and written to file. Then a line should follow with the two
              states for each transition. The ket and the bra orbitals are written
              followed by the singular values in the the usual format used by MOLCAS.
              </HELP>
              </KEYWORD>

:kword:`ORBItals`
  Print out the Molecular Orbitals read from each
  :file:`JOBIPH` file.

  .. xmldoc:: %%Keyword: Orbitals <advanced>
              Print out the molecular orbitals read from each JOBIPH file.

:kword:`OVERlaps`
  Print out the overlap integrals between the various orbital sets.

  .. xmldoc:: %%Keyword: Overlaps <advanced>
              Print out the overlap integrals between the various orbital sets.

:kword:`CIPRint`
  Print out the CI coefficients read from
  :file:`JOBIPH`.

  .. xmldoc:: %%Keyword: CIPrint <advanced>
              Print out the CI coefficients read from JOBIPH.

:kword:`THRS`
  The next line gives the threshold for printing CI coefficients. The
  default value is 0.05.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="THRS" APPEAR="Threshold for CI coeff." KIND="REAL" LEVEL="ADVANCED">
              %%Keyword: Thrs <advanced>
              <HELP>
              Enter the threshold for printing CI coefficients. Default 0.05.
              </HELP>
              </KEYWORD>

:kword:`DIPRrint`
  The next entry gives the threshold for printing dipole intensities.
  Default is 1.0D-5.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="DIPR" APPEAR="Threshold dipole intensities" KIND="REAL" MIN_VALUE="0.0" DEFAULT_VALUE="1.0D-5" LEVEL="ADVANCED">
              %%Keyword: DIPR <advanced>
              <HELP>
              Enter the threshold for printing dipole intensities.
              Default is 1.0D-5.
              </HELP>
              </KEYWORD>

:kword:`QIPRrint`
  The next entry gives the threshold for printing quadrupole intensities.
  Default is 1.0D-5.
  Will overwrite any value chosen for dipole intensities.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="QIPR" APPEAR="Threshold quadrupole intensities" KIND="REAL" MIN_VALUE="0.0" DEFAULT_VALUE="1.0D-5" LEVEL="ADVANCED">
              %%Keyword: QIPR <advanced>
              <HELP>
              Enter the threshold for printing quadrupole intensities.
              Default is 1.0D-5.
              Will overwrite any value choosen for dipole intensities.
              </HELP>
              </KEYWORD>

:kword:`QIALL`
  Print all quadrupole intensities.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="QIALL" APPEAR="Print all quadrupole intensities" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: QIALL <advanced>
              <HELP>
              Print all quadrupole intensities.
              </HELP>
              </KEYWORD>

:kword:`TMOS`
  Activate the computation of oscillators strengths (and transition moments) using the
  non-relativistic Hamiltonian with the explicit Coulomb-field vector operator (A) in
  the weak field approximation.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="TMOS" APPEAR="Transition moments - oscillator strength" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: TMOS <advanced>
              <HELP>
              Activate the computation of oscillators strengths (and transition moments) using the
              non-relativistic Hamiltonian with the explicit Coulomb-field vector operator (A) in
              the weak field approximation.
              </HELP>
              </KEYWORD>

:kword:`L-EFfective`
  Set the order of the Lebedev grids used in the interpolation of the solid angles
  in association with the :kword:`TMOS` option. Default value is 5.
  Other allowed values are: 7, 11, 17, 23, 29, 35, 41, 47, 53, and 59.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="L-EFF" APPEAR="Order of Lebedev integration" KIND="INT" MIN_VALUE="0" DEFAULT_VALUE="5" LEVEL="ADVANCED">
              %%Keyword: L-Eff <advanced>
              <HELP>
              Set the order of the Lebedev grids used in the interpolation of the solid angles
              in association with the TMOS option. Default value is 5.
              Other allowed values are: 7, 11, 17, 23, 29, 35, 41, 47, 53, and 59.
              </HELP>
              </KEYWORD>

:kword:`K-VEctor`
  Define the direction of the incident light for which we will
  compute transition moments and oscillator strengths. The keyword
  is followed by three reals specifying the direction. The values
  do not need to be normalized.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="K-VECTOR" APPEAR="The wave k-vector" KIND="REALS" SIZE="3" LEVEL="ADVANCED">
              %%Keyword: k-vector <advanced>
              <HELP>
              Define the direction of the incident light for which we will
              compute transition moments and oscillator strengths. The keyword
              is followed by three reals specifying the direction. The values
              do not need to be normalized.
              </HELP>
              </KEYWORD>

:kword:`RFPErt`
  :program:`RASSI` will read from :file:`RUNOLD` (if not present defaults to :file:`RUNFILE`) a response field contribution
  and add it to the Fock matrix.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="RFPE" APPEAR="Response field" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: Rfpe <basic>
              <HELP>
              RASSI will read from RUNOLD (if not present defaults to RUNFILE) a response field contribution
              and add it to the Fock matrix.
              </HELP>
              </KEYWORD>

:kword:`HCOM`
  The spin-free Hamiltonian is computed.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="HCOMPUTE" APPEAR="Computed Hamiltonian" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: HCom <basic>
              <HELP>
              The spin-free Hamiltonian is computed.
              </HELP>
              </KEYWORD>

:kword:`HEXT`
  It is read from the following few lines, as a triangular matrix: One element
  of the first row, two from the next, etc., as list-directed input of reals.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="HEXT" APPEAR="External Hamiltonian" KIND="STRINGS" LEVEL="ADVANCED">
              %%Keyword: HExt <advanced>
              <HELP>
              The spin-free Hamiltonian is read from a file instead of being computed.
              </HELP>
              It is read from the following entries, as a triangular matrix: One element
              of the first row, two from the next, etc., as list-directed input of reals.
              </KEYWORD>

:kword:`HEFF`
  A spin-free effective Hamiltonian is read from :file:`JOBIPH` instead of being computed.
  It must have been computed by an earlier program. Presently, this is done by
  a multi-state calculation using :program:`CASPT2`. In the future, other programs may add
  dynamic correlation estimates in a similar way. This keyword is not needed if the input
  file is in HDF5 format.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="HEFF" APPEAR="Effective Hamiltonian" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: HEff <advanced>
              <HELP>
              A spin-free effective Hamiltonian is read from JOBIPH instead of being computed.
              It must have been computed by an earlier program. Presently, this is done by
              a multi-state calculation using CASPT2.
              </HELP>
              </KEYWORD>

:kword:`EJOB`
  The spin-free effective Hamiltonian is assumed to be diagonal, with energies
  being read from a :file:`JOBIPH` or :file:`JOBMIX` file.
  If this keyword is used together with :kword:`HEFF`, or if the input file is
  an HDF5 file for which the effective Hamiltonian is automatically read, only
  the diagonal elements will be read and off-diagonal elements will be set to zero.
  This can be useful to use the SS-CASPT2 energies from a MS-CASTP2 calculation.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="EJOB" APPEAR="Read energies from file" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: EJob <advanced>
              <HELP>
              The spin-free effective Hamiltonian is assumed to be diagonal, with energies
              being read from a JOBIPH or JOBMIX file from e.g. a multi-state CASPT2 calculation.
              </HELP>
              </KEYWORD>

:kword:`TOFIle`
  Signals that a set of files with data from :program:`Rassi` should be
  created. This keyword is necessary if :program:`QmStat` is to be run
  afterwards.

  .. xmldoc:: %%Keyword: TOfile <basic>
              Signals that a set of files with data from Rassi should be
              created. This keyword is necessary if QmStat is to be run
              afterwards.

:kword:`XVIN`
  Demand for printing expectation values of all selected one-electron
  properties, for the input RASSCF wave functions.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="XVIN" APPEAR="input expectation values" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: XVIN <basic>
              <HELP>
              Demand printing expectation values of all selected one-electron
              properties, for the input RASSCF wave functions.
              </HELP>
              </KEYWORD>

:kword:`XVES`
  Demand for printing expectation values of all selected one-electron
  properties, for the spin-free eigenstates.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="XVES" APPEAR="spin-free expectation values" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: XVES <basic>
              <HELP>
              Demand printing expectation values of all selected one-electron
              properties, for the spin-free eigenstates.
              </HELP>
              </KEYWORD>

:kword:`XVSO`
  Demand for printing expectation values of all selected one-electron
  properties, for the spin--orbit states.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="XVSO" APPEAR="spin-orbit expectation values" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: XVSO <basic>
              <HELP>
              Demand printing expectation values of all selected one-electron
              properties, for the spin-orbit states.
              </HELP>
              </KEYWORD>

:kword:`EPRG`
  This computes the g matrix and principal g values for the
  states lying within the energy range supplied on the next line.
  A value of 0.0D0 or negative will select only the ground state,
  a value E will select all states within energy E of the ground state.
  The states should be ordered by increasing energy in the input.
  The angular momentum and spin--orbit coupling matrix elements
  need to be available (use keywords :kword:`SPIN` and :kword:`PROP`).
  For a more detailed description see ref :cite:`EPRG:2008`.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="EPRG" APPEAR="EPR g Matrix" KIND="REAL" LEVEL="ADVANCED" REQUIRE="SPIN,PROP">
              %%Keyword: EPRg <advanced>
              <HELP>
              This computes the g matrix and principal g values for the
              states lying within the energy range supplied on the next line.
              A value of 0.0D0 or negative will select only the ground state,
              a value E will select all states within energy E of the ground state.
              The states should be ordered by increasing energy in the input.
              The angular momentum and spin-orbit coupling matrix elements
              need to be available (use keywords SPIN and PROP).
              </HELP>
              </KEYWORD>

:kword:`MAGN`
  This computes the magnetic moment and magnetic susceptibility.
  On the next two lines you have to provide the magnetic field and
  temperature data. On the first line put the number of magnetic
  field steps, the starting field (in Tesla), size of the steps (in Tesla),
  and an angular resolution for sampling points in case of powder magnetization
  (for a value of 0.0d0 the powder magnetization is deactivated).
  The second line reads the number of temperature steps, the starting
  temperature (K), and the size of the temperature steps (K).
  The angular momentum and spin--orbit coupling matrix elements
  need to be available (use keywords :kword:`SPIN` and :kword:`PROP`).
  For a more detailed description see ref :cite:`MAGN:2009`.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="MAGN" APPEAR="Magnetism" KIND="REAL" LEVEL="DEV" REQUIRE="SPIN,PROP">
              %%Keyword: Magnetic properties <advanced>
              <HELP>
              This computes the magnetic moment and magnetic susceptibility.
              On the next two lines you have to provide the magnetic field and
              temperature data. On the first line put the number of magnetic
              field steps, the starting field (in Tesla), size of the steps (in Tesla),
              and an angular resolution for sampling points in case of powder magnetization
              (for a value of 0.0d0 the powder magnetization is deactivated).
              The second line reads the number of temperature steps, the starting
              temperature (K), and the size of the temperature steps (K).
              The angular momentum and spin-orbit coupling matrix elements
              need to be available (use keywords SPIN and PROP).
              For a more detailed description see J. Phys. Chem. A 113 6149.
              </HELP>
              </KEYWORD>

:kword:`HOP`
  Enables a trajectory surface hopping (TSH) algorithm which allow
  non-adiabatic transitions between electronic states during molecular
  dynamics simulation with :program:`DYNAMIX` program. The algorithm
  computes the scalar product of the amplitudes of different
  states in two consecutive steps. If the scalar product
  deviates from the given threshold a transition between the states
  is invoked by changing the root for the gradient computation.
  The current implementation is working only with SA-CASSCF.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="HOP" APPEAR="Trajectory surface hopping algorithm" KIND="SINGLE" LEVEL="DEV">
              %%Keyword: Hop <advanced>
              <HELP>
              Allows transitions between electronic states during molecular
              dynamics simulations.
              </HELP>
              </KEYWORD>

:kword:`STOVerlaps`
  Computes only the overlaps between the input states.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="STOV" APPEAR="State overlaps" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: StOverlaps <advanced>
              <HELP>
              Computes only the overlaps between the input states.
              </HELP>
              </KEYWORD>

:kword:`TRACk`
  Tries to follow a particular root during an optimization.
  Needs two :file:`JOBIPH` files (see :kword:`NrOfJobIphs`) with the same
  number of roots. The first file corresponds to the current iteration,
  the second file is the one from the previous iteration (taken as a reference).
  With this keyword :program:`RASSI` selects the root from the first :file:`JOBIPH`
  with highest overlap with the root that was selected in the previous
  iteration. It also needs :kword:`MDRlxRoot`, rather than :kword:`RlxRoot`,
  to be specified in :program:`RASSCF`.
  No other calculations are done by :program:`RASSI` when :kword:`Track`
  is specified.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="TRACK" APPEAR="Track root" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: Track <advanced>
              <HELP>
              Tries to follow a particular root during an optimization.
              </HELP>
              </KEYWORD>

:kword:`DQVD`
  Perfoms DQΦ diabatization :cite:`Hoyer:2014fk` by using properties that are computed with :program:`RASSI`.
  Seven properties must be computed with RASSI in order for this keyword to work
  (:math:`x`, :math:`y`, :math:`z`, :math:`xx`, :math:`yy`, :math:`zz`, :math:`1/r`), they will be automatically selected with the default input
  if the corresponding integrals are available (see keywords :kword:`MULT` and :kword:`EPOT` in :program:`GATEWAY`).
  At present, this keyword also requires :kword:`ALPHa` and :kword:`BETA`, where
  :kword:`ALPHa` is the parameter in front of :math:`rr` and :kword:`BETA` is the parameter
  in front of :math:`1/r`. When :kword:`ALPHa` and :kword:`BETA` are equal to zero, this
  method reduces to Boys localized diabatization :cite:`Subotnik:2008fk`.
  At present, this method only works for one choice of origin for each quantity.

  .. See Test/input/test393.input for an example.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="DQVD" APPEAR="DQV diabatization" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: DQVD <advanced>
              <HELP>
              Diabatizes by using dipole, quadrupole, and/or electrostatic potential.
              </HELP>
              </KEYWORD>

:kword:`ALPHa`
  :kword:`ALPHa` is the prefactor for the quadrupole term in DQΦ diabatization. This
  keyword must be used in conjunction with :kword:`DQVD` and :kword:`BETA`. You must
  specify a real number (e.g. :math:`\alpha = 1.0` not :math:`\alpha = 1`).

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="ALPHA" APPEAR="Alpha parameter" KIND="REAL" LEVEL="ADVANCED" REQUIRE="DQVD">
              %%Keyword: Alpha <advanced>
              <HELP>
              Alpha parameter in front of the quadrupole term for DQΦ diabatization.
              </HELP>
              </KEYWORD>

:kword:`BETA`
  :kword:`BETA` is the prefactor for the electrostatic potential term in DQΦ diabatization. This
  keyword must be used in conjunction with :kword:`DQVD` and :kword:`ALPHa`. You must
  specify a real number (e.g. :math:`\beta = 1.0` not :math:`\beta = 1`).

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="BETA" APPEAR="Beta parameter" KIND="REAL" LEVEL="ADVANCED" REQUIRE="DQVD">
              %%Keyword: Beta <advanced>
              <HELP>
              Beta parameter in front of the electrostatic potential term for DQΦ diabatization.
              </HELP>
              </KEYWORD>

:kword:`TRDI`
  Prints out the components and the module of the transition dipole
  vector. Only vectors with sizes large than 1.0D-4 a.u. are printed.
  See also the :kword:`TDMN` keyword.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="TRDI" APPEAR="Transition dipole" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: TRDI <advanced>
              <HELP>
              Prints out the components and the size of the transition dipole
              vector. Only vectors with sizes large than 1.0D-4 a.u. are printed.
              See also the TDMN keyword.
              </HELP>
              </KEYWORD>

:kword:`TDMN`
  Prints out the components and the module of the transition dipole
  vector. On the next line, the minimum size, in a.u., for the dipole
  vector to be printed must be given.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="TDMN" APPEAR="Transition dipole threshold" KIND="REAL" LEVEL="ADVANCED" REQUIRE="TRDI">
              %%Keyword: TDMN <advanced>
              <HELP>
              Prints out the components and the module of the transition dipole
              vector. On the next line, the minimum size, in a.u., for the dipole
              vector to be printed must be given.
              </HELP>
              </KEYWORD>

:kword:`TRD1`
  Prints the 1-electron (transition) densities to ASCII files and to
  the HDF5 file :file:`rassi.h5`.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="TRD1" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: TRD1 <advanced>
              <HELP>
              Prints the 1-electron (transition) densities to ASCII files and to
              the HDF5 file rassi.h5.
              </HELP>
              </KEYWORD>

:kword:`TRD2`
  Prints the 1/2-electron (transition) densities to ASCII files.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="TRD2" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: TRD2 <advanced>
              <HELP>
              Prints the 1/2-electron (transition) densities to ASCII files.
              </HELP>
              </KEYWORD>

:kword:`DYSOn`
  Enables calculation of Dyson amplitudes (an approximation of photo-electron intensities) between states that differ by exactly one in their number of electrons.

  Calculations are performed for spin-free states, and for spin-orbit coupled states if the keyword :kword:`SPINorbit` has also been specified. Note that spin-orbit coupled amplitudes are per default obtained from an approximation where a transformation is applied directly to the spin-free amplitudes rather than the Dyson orbitals, which may severly impact the accuracy. For a complete calculation also for spin-orbit states see the :kword:`DYSExport` keyword.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="DYSON" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: DYSON <advanced>
              <HELP>
              Enables calculation of Dyson amplitudes (an approximation of photo-electron intensities) between states that differ by exactly one in their number of electrons.
              </HELP>
              </KEYWORD>

:kword:`DYSExport`
  Requires the :kword:`DYSOn` keyword and enables exportation of Dyson orbitals (from which Dyson amplitudes are obtained). The next line specifies the number (starting from the first) of spin-free and spin-orbit states (two numbers, both mandatory) for which the exportation will be done. Note that the ordering of spin-free states depends on the ordering of JOBfiles, whereas spin-orbit states are always energy ordered.

  Dyson amplitudes for the spin-orbit states are here correctly obtained from a transformation of the Dyson orbitals (as opposed to the amplitudes, see :kword:`DYSOn` keywpord), but only for the specified number of initial states. Note that this calculation may be time consuming, i.e. the number of initial states should be limited.

  .. xmldoc:: <KEYWORD MODULE="RASSI" NAME="DYSEXPORT" KIND="INTS" SIZE="2" LEVEL="ADVANCED">
              %%Keyword: DYSEXPORT <advanced>
              <HELP>
              Requires the DYSOn keyword and enables exportation of Dyson orbitals (from which Dyson amplitudes are obtained). The next line specifies the number (starting from the first) of spin-free and spin-orbit states (two numbers, both mandatory) for which the exportation will be done. Note that the ordering of spin-free states depends on the ordering of JOBfiles, whereas spin-orbit states are always energy ordered.
              </HELP>
              </KEYWORD>

Input example
.............

::

  >>COPY  "Jobiph file 1" JOB001
  >>COPY  "Jobiph file 2" JOB002
  >>COPY  "Jobiph file 3" JOB003

  &RASSI
  NR OF JOBIPHS= 3 4 2 2    --- 3 JOBIPHs. Nr of states from each.
  1 2 3 4; 3 4; 3 4         --- Which roots from each JOBIPH.
  CIPR; THRS= 0.02
  Properties= 4; 'MltPl  1'  1   'MltPl  1'  3    'Velocity'  1 'Velocity'  3
  * This input will compute eigenstates in the space
  * spanned by the 8 input functions. Assume only the first
  * 4 are of interest, and we want natural orbitals out
  NATO= 4

.. xmldoc:: </MODULE>
