.. index::
   single: Program; MBPT2
   single: MBPT2

.. _UG\:sec\:mbpt2:

:program:`mbpt2`
================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="MBPT2">
            %%Description:
            <HELP>
            This program computes the second order Many Body Perturbation Theory
            correction to an SCF wavefunction.
            </HELP>

.. _UG\:sec\:mbpt2_description:

Description
-----------

The :program:`MBPT2` program of the |molcas| program system computes
the second order correlation energy and the reference weight for a
closed-shell Hartree--Fock reference wave function, based on a
Møller--Plesset partitioning of the Hamiltonian and canonical orbitals.

If :program:`SEWARD` performed a Cholesky decomposition of the two-electron integrals prior to running
the :program:`MBPT2` program, Cholesky vectors will be employed for computing
the second order energy correction. This is done by first transforming the
Cholesky vectors to MO basis and subsequently computing the :math:`(ai|bj)` integrals.
These integrals are either computed, stored on disk, and then read back in to
memory during the energy evaluation (i.e. mimicking a conventional calculation)
or they may be computed on-the-fly. The user may choose either algorithm
through the Cholesky-specific options described below.

If :program:`SEWARD` did not perform a Cholesky decomposition,
the transformation of the two-electron integrals in AO basis
(:math:`\mathcal{O}(N^4)`, where :math:`N` is the number of basis functions)
to the exchange operator matrices :math:`\mat{K}^{ij}` in MO basis
(:math:`\mathcal{O}(O^2)` matrices of size :math:`V^2`, where :math:`O` and :math:`V` denote the number
of occupied and virtual orbitals, respectively), is either done
conventionally, using the two-electron integral file :file:`ORDINT`, which
was generated in a previous step by the :program:`SEWARD` integral code.

.. _UG\:sec\:mbpt2_dependencies:

Dependencies
------------

The :program:`MBPT2` program requires the communications file :file:`RUNFILE`.
It contains specifications processed by :program:`SEWARD`,
the Hartree--Fock canonical orbitals, eigenvalues and energy generated
by :program:`SCF`.
For Cholesky-based calculations, all Cholesky related files (see the
manual pages for :program:`SEWARD`) must be available, whereas
for conventional (not integral-direct) calculations
the two-electron integral file :file:`ORDINT`
is required. Hence, before running :program:`MBPT2`, a :program:`SEWARD`
and a :program:`SCF` run have to be performed.

.. index::
   pair: Files; MBPT2

.. _UG\:sec\:mbpt2_files:

Files
-----

Input files
...........

:program:`MBPT2` will use the following input
files: :file:`ONEINT`, :file:`ORDINT`, :file:`RUNFILE`.
For Cholesky runs: :file:`CHVEC`, :file:`CHORST`, :file:`CHRED` and
:file:`CHOR2F`
(for more information see :numref:`UG:sec:files_list`).

.. Intermediate files
   ..................

   All the intermediate files are created, used and removed
   automatically, unless you yourself create a link or a file
   with the specified name.

   .. class:: filelist

   :file:`MOLINT*`
     Resulting file of transformed integrals.
     Scratch file; conventional calculation only.

   :file:`LUHLFn*`
     :math:`n`\=1 to 3. Intermediate files used in the 1st, 2nd, and 3rd, respectively,
     transformation step. Conventional calculation only.

Output files
............

.. class:: filelist

:file:`RUNFILE`
  File for communication of auxiliary information.

.. _UG\:sec\:mbpt2_input:

Input
-----

Below follows a description of the input to :program:`MBPT2`.
The input for each module is preceded by its name like: ::

  &MBPT2

No compulsory keywords are required for :program:`MBPT2`.
The reference statement mentioned
above is sufficient for a default :program:`MBPT2` run.

Optional keywords
.................

.. class:: keywordlist

:kword:`TITLe`
  The line following this line is regarded as a title line

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="TITLE" KIND="STRING" LEVEL="BASIC">
              <HELP>
              Print a title line
              </HELP>
              %%Keyword: Title <basic>
              The line following this line is regarded as a title line
              </KEYWORD>

:kword:`PRINt`
  Specifies the general print level of the calculation. An integer
  has to be supplied as argument.
  The default value, 0, is recommended for production calculations.

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="PRINT" APPEAR="Print level" LEVEL="ADVANCED" KIND="INT" DEFAULT_VALUE="0">
              %%Keyword: Print <advanced>
              <HELP>
              Specifies the general print level of the calculation. An integer
              has to be supplied as argument.
              </HELP>
              The default value, 0, is recommended for production calculations.
              </KEYWORD>

:kword:`FREEze`
  Specifies the total number of frozen occupied orbitals.
  The lowest-energy occupied orbitals are then automatically identified and frozen.
  The keyword takes as argument one integer.
  Incompatible with the :kword:`FROZen` keyword.

  .. xmldoc:: <SELECT MODULE="MBPT2" NAME="ORBITAL_FREEZE" APPEAR="Frozen orbitals selection" CONTAINS="FREEZE,FROZEN">

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="FREEZE" APPEAR="freeze orbitals" LEVEL="ADVANCED" KIND="INT" DEFAULT_VALUE="0" EXCLUSIVE="FROZEN">
              %%Keyword: Freeze <advanced>
              <HELP>
              Specifies the total number of frozen occupied orbitals.
              </HELP>
              Incompatible with the FROZen keyword.
              </KEYWORD>

:kword:`FROZen`
  Specifies the number of frozen occupied orbitals in each of the irreducible
  representations (irreps) of the subgroup of :math:`D_{2h}` in which the system
  is represented. The counting of the orbitals follows the *increasing*
  orbital energy within each irrep, with those orbitals being frozen first
  that correspond to lowest orbital energies.
  The keyword takes as argument *nIrrep* (# of irreps) integers.
  Incompatible with the :kword:`FREEze` keyword.
  Default is to freeze non-valence orbitals.

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="FROZEN" APPEAR="frozen orbitals" LEVEL="ADVANCED" KIND="INTS_LOOKUP" SIZE="NSYM" EXCLUSIVE="FREEZE">
              %%Keyword: Frozen <advanced>
              <HELP>
              Specifies the number of frozen orbitals in each irrep
              of the point group. The orbitals with the lowest energies are frozen.
              </HELP>
              The keyword takes as argument nIrrep (# of irreps) integers.
              Incompatible with the FREEze keyword.
              Default is to freeze non-valence orbitals.
              </KEYWORD>

  .. xmldoc:: </SELECT>

:kword:`DELEted`
  Specifies the number of deleted orbitals in each of the irreducible
  representations (irreps) of the subgroup of :math:`D_{2h}` in which the system
  is represented. The counting of the orbitals follows the *decreasing*
  orbital energy within each irrep, with those orbitals being deleted first
  that correspond to highest orbital energies.
  The keyword takes as argument *nIrrep* (# of irreps) integers.

  **NOTE:** Those orbitals, which have been deleted already in the
  :program:`SCF` calculation (cf. :kword:`SPDElete`, :kword:`OVLDelete` of
  the :program:`SCF` program description) are never seen by the
  :program:`MBPT2` program and hence are **not** to be deleted again with
  the present option.

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="DELETED" APPEAR="deleted orbitals" LEVEL="ADVANCED" KIND="INTS_LOOKUP" SIZE="NSYM">
              %%Keyword: Deleted <advanced>
              <HELP>
              Specifies the number of deleted orbitals in each irrep of the point group.
              The orbitals with the highest energies are deleted.
              </HELP>
              The keyword takes as argument nIrrep (# of irreps) integers.
              NOTE: Those orbitals, which have been deleted already in the
              SCF calculation (cf. SPDElete, OVLDelete of
              the SCF program description) are never seen by the
              MBPT2 program and hence are not to be deleted again with
              the present option.
              </KEYWORD>

:kword:`SFROzen`
  Allows to specify specific orbitals to freeze in each of the irreducible
  representations (irreps) of the subgroup of :math:`D_{2h}` in which the system
  is represented. In the 1st line after the keyword the number of orbitals
  to freeze for each irrep is specified (*nIrrep* (# of irreps) integers).
  The next :math:`\leq` *nIrrep* lines reference the orbitals to freeze for the
  related irrep, following an enumeration of the individual orbitals
  of 1, 2, 3,... according to
  *increasing* orbital energy. Note that the orbital reference numbers
  obey the original ordering and also include those orbitals which
  may have been frozen already by the
  :kword:`FROZen` or :kword:`FREEze` options. If the corresponding irrep does not contain any
  specific orbitals to freeze (i.e. a zero was supplied for this irrep in the
  1st line), no line orbital reference input line is supplied for that irrep.

  .. xmldoc:: %%Keyword: Sfrozen <advanced>
              Allows to specify specific orbitals to freeze in each of the irreducible
              representations (irreps) of the subgroup of D2h in which the system
              is represented. In the 1st line after the keyword the number of orbitals
              to freeze for each irrep is specified (nIrrep (# of irreps) integers).
              The next <= nIrrep lines reference the orbitals to freeze for the
              related irrep, following an enumeration of the individual orbitals
              of 1, 2, 3,... according to
              increasing orbital energy. Note that the orbital reference numbers
              obey the original ordering and also include those orbitals which
              may have been frozen already by the
              FROZEN option. If the corresponding irrep does not contain any
              specific orbitals to freeze (i.e. a zero was supplied for this irrep in the
              1st line), no line orbital reference input line is supplied for that irrep.

:kword:`SDELeted`
  Allows to specify specific orbitals to delete in each of the irreducible
  representations (irreps) of the subgroup of :math:`D_{2h}` in which the system
  is represented. In the 1st line after the keyword the number of orbitals
  to delete for each irrep is specified (*nIrrep* (# of irreps) integers).
  The next :math:`\leq` *nIrrep* lines reference the orbitals to delete for the
  related irrep, following an enumeration of the individual orbitals
  of 1, 2, 3,... according to
  *increasing* orbital energy. Note that the orbital reference numbers
  obey the original ordering.
  If the corresponding irrep does not contain any
  specific orbitals to freeze (i.e. a zero was supplied for this irrep in the
  1st line), no line orbital reference input line is supplied for that irrep.

  .. xmldoc:: %%Keyword: Sdeleted <advanced>
              Allows to specify specific orbitals to delete in each of the irreducible
              representations (irreps) of the subgroup of D2h in which the system
              is represented. In the 1st line after the keyword the number of orbitals
              to delete for each irrep is specified (nIrrep (# of irreps) integers).
              The next <= nIrrep lines reference the orbitals to delete for the
              related irrep, following an enumeration of the individual orbitals
              of 1, 2, 3,... according to
              increasing orbital energy. Note that the orbital reference numbers
              obey the original ordering.
              If the corresponding irrep does not contain any
              specific orbitals to freeze (i.e. a zero was supplied for this irrep in the
              1st line), no line orbital reference input line is supplied for that irrep.

:kword:`GHOStdelete`
  Excludes from PT2 treatment orbitals localized on ghost atoms. A threshold for this selection must be specified.

  .. xmldoc:: %%Keyword: GHOS <advanced>
              Excludes from PT2 treatment orbitals localized on ghost atoms. A threshold for this selection must be specified.

:kword:`LUMOrb`
  Molecular orbital coefficients and energies read from :file:`INPORB` file rather
  than :file:`RunFile`.

  .. xmldoc:: %%Keyword: LUMO <basic>
              Molecular orbital coefficients and energies read from INPORB file rather
              than RunFile.

:kword:`EREF`
  Specifies the value of the reference energy. Available only in combination
  with :kword:`LumOrb`. Default value of the reference energy is set to zero.

  .. xmldoc:: %%Keyword: EREF <basic>
              Specifies the value of the reference energy. Available only in combination
              with LumOrb. Default value of the reference energy is set to zero.

:kword:`TEST`
  If this keyword is specified the input is checked without performing any
  calculation.

  .. xmldoc:: %%Keyword: TEST <basic>
              If this keyword is specified the input is checked without performing any
              calculation.

:kword:`T1AM`
  Singles amplitudes/energy introduced according to Thouless formula.
  An INPORB file containing MOs different from HF orbitals is required.

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="T1AM" APPEAR="Thouless singles amplitudes" LEVEL="ADVANCED" KIND="SINGLE">
              %%Keyword: T1AM <advanced>
              <HELP>
              Singles amplitudes/energy introduced according to Thouless formula.
              An INPORB file containing MOs different from HF orbitals is required.
              </HELP>
              </KEYWORD>

:kword:`LOVMp2`
  "Freeze-and-Delete" type of MP2, available only in connection with Cholesky or RI.
  An example of input for the keyword :kword:`LOVM` is the following: ::

    LovMP2
    2  0.2  (nCenters,thrs)
    C1 N    (Center labels)
    DoMP2

  In this case, both occupied and virtual orbitals (localized by the program) are divided in two groups: those (A) mainly located on the
  two (symmetry independent) centers C1 and C2, and the remaining ones (B), which are obviously "outside" this region.
  The value of the threshold (between 0 and 1) is used to perform this selection
  (in the example, 20% of the gross Mulliken population of a given orbital on the specified atoms).
  By default, the MP2 calculation is performed only for the correlating orbitals associated with the region A ("active site").
  The keyword :kword:`DoMP2` is optional and forces the program to perform also an independent MP2 calculation on
  the "frozen region" (B).
  Alternatively, one can specify the keyword :kword:`VirAll` in order to use all virtual orbitals as correlating space for the
  occupied orbitals of the active site.

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="LOVMP2" APPEAR="Localized occupied-virtual MP2" LEVEL="ADVANCED" KIND="CUSTOM">
              %%Keyword: LOVM <advanced>
              <HELP>
              "Freeze-and-Delete" type of MP2, available only in connection with Cholesky or RI.
              An example of input for the keyword LOVM is the following:
              ||
              ||LovMP2
              ||2  0.2  (nCenters,thrs)
              ||C1 N    (Center labels)
              ||DoMP2
              ||
              In this case, both occupied and virtual orbitals (localized by the program) are divided in two groups: those (A) mainly located on the
              two (symmetry independent) centers C1 and N, and the remaining ones (B), which are obviously "outside" this region.
              The value of the threshold (between 0 and 1) is used to perform this selection
              (in the example, 20% of the gross Mulliken population of a given orbital on the specified atoms).
              By default, the MP2 calculation is performed only for the correlating orbitals associated with the region A ("active site").
              The keyword DoMP2 is optional and forces the program to perform also an independent MP2 calculation on
              the "frozen region" (B).
              Alternatively, one can specify the keyword VirAll in order to use all virtual orbitals as correlating space for the
              occupied orbitals of the active site.
              </HELP>
              </KEYWORD>

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="DOMP2" LEVEL="UNDOCUMENTED" KIND="SINGLE" />

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="VIRALL" LEVEL="UNDOCUMENTED" KIND="SINGLE" />

:kword:`FNOMp2`
  Performs a Frozen Natural Orbital (FNO) MP2 calculation, available only in combination with Cholesky or RI integral representation.
  An example of input for the keyword :kword:`FNOM` is the following: ::

    FNOMp2
    0.4
    DoMP2

  The keyword :kword:`FNOM` has one compulsory argument (real number in ]0,1]) specifying the fraction of virtual orbitals
  (in each irrep) to be retained in the FNO-MP2 calculation.
  The keyword :kword:`DoMP2` is optional and used to compute the (estimated) correction for the truncation error.

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="FNOMP2" APPEAR="Frozen natural orbital MP2" LEVEL="ADVANCED" KIND="REAL">
              %%Keyword: FNOM <advanced>
              <HELP>
              Performs a Frozen Natural Orbital (FNO) MP2 calculation, available only in combination with Cholesky or RI integral representation
              An example of input for the keyword FNOM is the following:
              ||
              ||FNOMp2
              || 0.4
              ||DoMP2
              ||
              The keyword FNOM has one compulsory argument (real number in ]0,1]) specifying the fraction of virtual orbitals
              (in each irrep) to be retained in the FNO-MP2 calculation.
              The keyword DoMP2 is optional and used to compute the (estimated) correction for the truncation error.
              </HELP>
              </KEYWORD>

:kword:`PRPT`
  Multipole moments (dipoles and quadrupoles) are calculated and printed. The moments
  are calculated by using a variational one-particle MP2 density matrix.
  The calculation of the density matrix substantially increases
  the computational effort compared to an ordinary energy calculation. If the call
  to :program:`MBPT2` is followed by a :program:`LoProp` call the variational MP2
  density matrix will automatically be passed on to that module when this keyword
  is active.

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="PRPT" APPEAR="Properties" LEVEL="BASIC" KIND="SINGLE">
              <HELP>
              Multipole moments are calculated and printed.
              </HELP>
              %%Keyword: PrPt <basic>
              Multipole moments (dipoles and quadrupoles) are calculated and printed.
              The computational effort is increased substantially compared to an energy-only
              calculation.
              </KEYWORD>

:kword:`GRDT`
  Variational one and two-particle MP2 densities are calculated to prepare for
  analytical gradient calculations. The default for subsequent gradient
  calculations are changed from numerical to analytical when this keyword is
  invoked. When using :program:`mbpt2` in a :program:`slapaf`\-loop with only :math:`C_1` symmetry
  analytical gradients are automatically default and this keyword is not
  needed. :kword:`grdt`
  prints Multipole moments and prepare for :program:`LoProp` in the exact same way
  as :kword:`prpt`.
  Use of this keyword therefore makes it
  redundant (but harmless) to also specify the keyword :kword:`prpt`.

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="GRDT" APPEAR="Analytic Gradient" LEVEL="BASIC" KIND="SINGLE" EXCLUSIVE="NOGRDT">
              %%Keyword: Grdt <advanced>
              <HELP>
              Analytical gradients are used in subsequent gradient calculations.
              </HELP>
              </KEYWORD>

:kword:`NOGRdt`
  Disables the calculation of variational densities for analytical gradients.
  This is useful to cancel the implicit :kword:`grdt` added when using :program:`mbpt2`
  inside a :program:`slapaf`\-loop, if no analytical gradients are actually needed.
  Note that using the :kword:`Numerical` keyword in :program:`gateway` already disables
  :kword:`grdt`, so :kword:`nogrdt` is only needed in some advanced situations.

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="NOGRDT" APPEAR="No Analytic Gradient" LEVEL="BASIC" KIND="SINGLE" EXCLUSIVE="GRDT">
              %%Keyword: NoGrdt <advanced>
              <HELP>
              Disables calculation of variational densities for analytical gradients.
              </HELP>
              </KEYWORD>

Optional keywords specific to Cholesky calculations
...................................................

*Observe* that these keywords are disregarded if the integrals
were not Cholesky decomposed by :program:`SEWARD`. Furthermore, they
are disregarded for algorithm 0 (see below).

.. class:: keywordlist

:kword:`CHOAlgorithm`
  Takes as argument one positive integer specifying
  the algorithm to use for Cholesky MP2.
  Options: 0 [generate MO integrals on disk from Cholesky vectors],
  1 [compute integrals on-the-fly, minimal operation count, level 2 BLAS],
  2 [compute integrals on-the-fly, not minimal operation count, level 3 BLAS],
  Default is 2.

  .. xmldoc:: <GROUP MODULE="MBPT2" NAME="CHOINPUT" APPEAR="Cholesky input section" KIND="BOX" LEVEL="ADVANCED">

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="CHOALGORITHM" APPEAR="CD algorithm" LEVEL="ADVANCED" KIND="INT" DEFAULT_VALUE="2">
              %%Keyword: ChoAlgorithm <advanced>
              <HELP>
              Specifies the algorithm to use for Cholesky MP2.
              </HELP>
              Options:
              ||0 [generate MO integrals on disk from Cholesky vectors]
              ||1 [compute integrals on-the-fly, minimal operation count]
              ||2 [compute integrals on-the-fly, minimal disk access (default)]
              </KEYWORD>

:kword:`VERBose`
  Increases printing from the Cholesky MP2 routines, although not
  by much.
  Default is (almost) no printing.

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="VERBOSE" APPEAR="Verbose printout" LEVEL="ADVANCED" KIND="SINGLE">
              %%Keyword: Verbose <advanced>
              <HELP>
              Increases printing from the Cholesky MP2 routines.
              </HELP>
              </KEYWORD>

:kword:`DECOmpose`
  Requests Cholesky decomposition of the :math:`(ai|bj)` integrals.
  Unless user-defined (see below), the threshold used is identical
  to that used by :program:`SEWARD` for decomposing the two-electron
  integrals. Default is to not decompose.

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="DECOMPOSE" APPEAR="MO integrals CD" LEVEL="ADVANCED" KIND="SINGLE">
              %%Keyword: Decompose <advanced>
              <HELP>
              Requests Cholesky decomposition of the (ai|bj) integrals.
              </HELP>
              </KEYWORD>

:kword:`THRCholesky`
  Specifies the threshold for :math:`(ai|bj)` Cholesky decomposition.
  Takes as argument one real number.
  Default is the threshold used by :program:`SEWARD` for decomposing the two-electron
  integrals.

  .. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="THRCHOLESKY" APPEAR="CD threshold" LEVEL="ADVANCED" KIND="REAL">
              %%Keyword: ThrCholesky <advanced>
              <HELP>
              Specifies the threshold for (ai|bj) Cholesky decomposition.
              </HELP>
              </KEYWORD>

  .. xmldoc:: </GROUP>

:kword:`NODEcompose`
  Turns off Cholesky decomposition of the :math:`(ai|bj)` integrals.
  Default is to not decompose.

  .. xmldoc:: %%Keyword: Nodecompose <advanced>
              Turns off Cholesky decomposition of the (ai|bj) integrals.

:kword:`SPAN`
  Specifies the span factor used for :math:`(ai|bj)` Cholesky decomposition.
  Takes as argument one real number.
  Default is the span factor used by :program:`SEWARD` for decomposing the two-electron
  integrals.

  .. xmldoc:: %%Keyword: Span <advanced>
              Specifies the span factor used for (ai|bj) Cholesky decomposition.

:kword:`MXQUal`
  Specifies the max. number of qualified diagonals treated during :math:`(ai|bj)` Cholesky decomposition.
  Takes as argument one integer.
  Default is 10% of the max. rank of :math:`(ai|bj)`, although never more than 200.

  .. xmldoc:: %%Keyword: MxQual <advanced>
              Specifies the max. number of qualified diagonals treated during (ai|bj) Cholesky decomposition.

:kword:`PRESort`
  Presort the MO Cholesky vectors according to the batches over occupied orbitals.
  This will reduce the amount of I/O performed during on-the-fly
  assembly of the :math:`(ai|bj)` integrals.
  This keyword is obsolete.

  .. Default is to sort when more than 2 batches over occupied orbitals are required.

  .. xmldoc:: %%Keyword: Presort <advanced>
              Presort the MO Cholesky vectors according to the batches over occupied orbitals.

Limitations
...........

The maximum number of selectively frozen
:kword:`SFRO` or selectively deleted orbitals
:kword:`SDEL` in each symmetry is limited to 50.

The limitations on the number of basis functions are the same as specified
for :program:`SEWARD`.

Input example
.............

::

  &MBPT2
  Title
   H2O:  O(9.5/4.2), H(4/2)
  * The lowest energy occupied orbital in the repr. no.1 will be frozen in
  * MBPT2 calculations. The number of representations is 4 and all zeros
  * must be explicitly given
  Frozen
  1 0 0 0
  * Two highest energy external orbitals in the repr. no.3 will be deleted
  * in MBPT2 calculations. The number of representations is 4 and all
  * zeros must be explicitly given
  Deleted
  0 0 2 0
  * One occupied orbital in symmetry no.1 will be additionally frozen by
  * using the SFRO option. Let it be the third SCF occupied orbital in
  * this symmetry
  sFrozen
  1 0 0 0   // Gives the number of frozen orbitals in each symmetry
  3         // Gives the frozen orbital reference number in symmetry no. 1
  * Two external orbitals in symmetry no.1 and one external orbital in
  * symmetry 3 will be deleted. In symmetry 1 let it be the second and
  * third external orbitals, and in symmetry 3 the third (already deleted
  * in by using the option DELE) external orbital
  sDeleted
  2 0 1 0   // Gives the number of orbitals to be deleted in each symmetry
  2 3       // Gives the reference numbers of external orbitals in sym. 1
  3         // Gives the reference number of the external orb. in sym. 3

.. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="SOSMP2" LEVEL="UNDOCUMENTED" KIND="SINGLE" />

.. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="FORCEBATCH" LEVEL="UNDOCUMENTED" KIND="SINGLE" />

.. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="OEDTHRESHOLD" LEVEL="UNDOCUMENTED" KIND="REAL" />

.. xmldoc:: <KEYWORD MODULE="MBPT2" NAME="OSFACTOR" LEVEL="UNDOCUMENTED" KIND="REAL" />

.. xmldoc:: </MODULE>
