.. index::
   single: Program; LOCALISATION
   single: LOCALISATION

.. _UG\:sec\:localisation:

:program:`localisation`
=======================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. _UG\:sec\:localisation_description:

Description
-----------

.. xmldoc:: <MODULE NAME="LOCALISATION">
            %%Description:
            <HELP>
            The LOCALISATION program of the molcas program system generates
            localised occupied orbitals according to one of the following procedures:
            Pipek-Mezey, Boys, Edmiston-Ruedenberg, or Cholesky.
            Orthonormal, linearly independent, local virtual orbitals may also be
            generated from projected atomic orbitals (Cholesky PAOs).
            </HELP>

The :program:`LOCALISATION` program of the |molcas| program system generates
localised occupied orbitals according to one of the following procedures:
Pipek--Mezey :cite:`Pipek:89`,
Boys :cite:`Boys:60,Foster:60`,
Edmiston--Ruedenberg :cite:`Edmiston:63`, or
Cholesky :cite:`Aquilante:06a`.
Orthonormal, linearly independent, local orbitals may also be
generated from projected atomic orbitals (Cholesky PAOs) :cite:`Aquilante:06a`.

.. compound::

  Orbital localisation makes use of the fact that a Hartree-Fock wave function
  is invariant under unitary transformations of the occupied orbitals,

  .. math:: \tilde{C}_{\mu i} = \sum_j C_{\mu j} \mat{U}_{ji} ,

  where :math:`\mat{U}` is unitary (i.e. orthogonal for real orbitals).
  The same is true for the inactive or active orbitals in a CASSCF wave function.
  Whereas the Pipek--Mezey :cite:`Pipek:89`,
  Boys :cite:`Boys:60,Foster:60`, and
  Edmiston--Ruedenberg :cite:`Edmiston:63` procedures define :math:`\mat{U}`
  through an iterative maximization of a localisation functional,
  the Cholesky orbitals are simply defined through the Cholesky decomposition
  of the one-electron density, i.e.

  .. math:: \sum_i \tilde{C}_{\mu i}\tilde{C}_{\nu i} = P_{\mu\nu} = \sum_i C_{\mu i} C_{\mu i} .

  Cholesky orbitals are thus not optimum localised orbitals by any of the
  Pipek--Mezey, Boys, or Edmiston--Ruedenberg measures, but rather inherit locality
  from the density matrix, see :cite:`Aquilante:06a` for details.

Although these localisation schemes are mostly meant for localising occupied
orbitals (except for PAOs which are defined for the virtual orbitals), the
:program:`LOCALISATION` program will attempt to localise any set of orbitals
that the user specifies. This means that it is possible to mix
occupied and virtual orbitals and thereby break the Hartree--Fock
invariance. The default settings, however, do not break the invariance.

For Pipek--Mezey, Boys, and Edmiston--Ruedenberg localisations, iterative
optimizations are carried out. We use
the :math:`\eta`-steps of Subotnik *et al.* :cite:`Subotnik:04` for
Edmiston--Ruedenberg, whereas the traditional Jacobi sweeps (consecutive
two-by-two orbital rotations) :cite:`Pipek:89,Subotnik:04`
are employed for the Pipek--Mezey and Boys schemes.

.. _UG\:sec\:localisation_dependencies:

Dependencies
------------

The :program:`LOCALISATION` program requires the one-electron integral file
:file:`ONEINT` and the communications file :file:`RUNFILE`,
which contains, among other data, the
basis set specifications processed by :program:`GATEWAY` and :program:`SEWARD`.
In addition, the Edmiston--Ruedenberg procedure requires the presence
of Cholesky decomposed two-electron integrals produced by :program:`SEWARD`.

.. index::
   pair: Files; LOCALISATION

.. _UG\sec\:localisation_files:

Files
-----

Below is a list of the files that are used/created by the program
:program:`LOCALISATION`.

Input files
...........

:program:`LOCALISATION` will use the following input
files: :file:`ONEINT`, :file:`RUNFILE`, :file:`INPORB`.
For Edmiston--Ruedenberg localisation,
it also needs :file:`CHVEC`, :file:`CHRED` and :file:`CHORST` files
(for more information see :numref:`UG:sec:files_list`).

Output files
............

.. class:: filelist

:file:`LOCORB`
  Localised orthonormal orbital output file.
  Note that :file:`LOCORB` contains all orbitals (localised as well as non-localised
  according to the input specification).

:file:`DPAORB`
  Linearly dependent nonorthonormal projected atomic orbital output file
  (only produced for PAO runs).

:file:`IPAORB`
  Linearly independent nonorthonormal projected atomic orbital output file
  (only produced for PAO runs).

:file:`RUNFILE`
  Communication file for subsequent programs.

:file:`MD_LOC`
  Molden input file for molecular orbital analysis.

.. index::
   pair: Input; LOCALISATION

.. _UG\:sec\:localisation_input:

Input
-----

Below follows a description of the input to :program:`LOCALISATION`.
The :program:`LOCALISATION` program section of the |molcas| input is bracketed by
a preceding program reference ::

  &LOCALISATION

Optional general keywords
.........................

.. class:: keywordlist

:kword:`FILEorb`
  The next line specifies the filename containing the input orbitals that will
  be localised. By default a file named :file:`INPORB` will be used.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="FILE" APPEAR="Orbitals file" KIND="STRING" LEVEL="BASIC">
              %%Keyword: FileOrb <basic>
              <HELP>
              The next line specifies the filename containing the input orbitals that will
              be localised. By default a file named INPORB will be used.
              </HELP>
              </KEYWORD>

:kword:`NORBitals`
  The following line specifies the number of orbitals to localise in each
  irreducible representation. The default is to localise all occupied
  orbitals as specified in the :file:`INPORB` input file, except for PAO runs where
  all the virtual orbitals are treated by default.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="NORB" APPEAR="Number of orbitals" LEVEL="BASIC" KIND="INTS_LOOKUP" SIZE="NSYM">
              <HELP>
              Please, specify the number of orbitals to localise in each irrep.
              </HELP>
              %%Keyword: NORB <basic>
              The following line specifies the number of orbitals to localise in each
              irreducible representation. The default is to localise all occupied
              orbitals as specified in the INPORB input file, except for PAO runs where
              all the virtual orbitals are treated by default.
              </KEYWORD>

:kword:`NFROzen`
  The following line specifies the number of orbitals to freeze in each
  irreducible representation. The default is not to freeze any orbitals,
  except for the localisations of the virtual space (see keywords :kword:`PAO` and
  :kword:`VIRTual`) where the default is to freeze all occupied orbitals (occupation
  number different from zero, as reported in the :file:`INPORB` file).

  .. xmldoc:: <SELECT MODULE="LOCALISATION" NAME="ORBITAL_FREEZE" APPEAR="Frozen orbitals selection" CONTAINS="NFROZEN,FREEZE">

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="NFROZEN" APPEAR="Orbitals to freeze" LEVEL="BASIC" KIND="INTS_LOOKUP" SIZE="NSYM" EXCLUSIVE="FREEZE">
              <HELP>
              Please, specify the number of orbitals to freeze in each irrep.
              </HELP>
              %%Keyword: NFRO <basic>
              The following line specifies the number of orbitals to freeze in each
              irreducible representation. The default is not to freeze any orbitals,
              except for the localisations of the virtual space (see keywords PAO and
              VIRTual) where the default is to freeze all occupied orbitals (occupation
              number different from zero, as reported in the INPORB file).
              </KEYWORD>

:kword:`FREEze`
  Implicit frozen core option. The default is not to freeze any orbitals,
  except for the localisations of the virtual space (see keywords :kword:`PAO` and
  :kword:`VIRTual`) where the default is to freeze all occupied orbitals (occupation
  number different from zero, as reported in the :file:`INPORB` file).
  The definition of core orbitals is taken from program :program:`SEWARD`.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="FREEZE" APPEAR="Freeze core orbitals" LEVEL="BASIC" KIND="SINGLE" EXCLUSIVE="NFROZEN">
              <HELP>
              Freeze the core orbitals as defined by SEWARD.
              </HELP>
              %%Keyword: FREE <basic>
              Implicit frozen core option. The default is not to freeze any orbitals,
              except for the localisations of the virtual space (see keywords PAO and
              VIRTual) where the default is to freeze all occupied orbitals (occupation
              number different from zero, as reported in the INPORB file).
              </KEYWORD>

  .. xmldoc:: </SELECT>

:kword:`OCCUpied`
  Requests that the occupied orbitals should be localised. This is the default
  except for PAO where the default is virtual.

  .. xmldoc:: <SELECT MODULE="LOCALISATION" NAME="LOC_ORB" APPEAR="Orbitals to localise" CONTAINS="OCCU,VIRT,ALL">

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="OCCU" APPEAR="Localise occupied orbitals" LEVEL="BASIC" KIND="SINGLE" EXCLUSIVE="VIRT,ALL">
              %%Keyword: OCCU <basic>
              <HELP>
              Requests that the occupied orbitals should be localised.
              </HELP>
              This is the default except for PAO where the default is virtual.
              </KEYWORD>

:kword:`VIRTual`
  Requests that the virtual orbitals should be localised. The default is
  to localise the occupied orbitals, except for PAO where the default is
  virtual.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="VIRT" APPEAR="Localise virtual orbitals" LEVEL="BASIC" KIND="SINGLE" EXCLUSIVE="OCCU,ALL">
              %%Keyword: VIRT <basic>
              <HELP>
              Requests that the virtual orbitals should be localised.
              </HELP>
              The default is
              to localise the occupied orbitals, except for PAO where the default is
              virtual.
              </KEYWORD>

:kword:`ALL`
  Requests that all orbitals should be localised. The default is
  to localise the occupied orbitals, except for PAO where the default is
  virtual.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="ALL" APPEAR="Localise all orbitals" LEVEL="BASIC" KIND="SINGLE" EXCLUSIVE="OCCU,VIRT">
              %%Keyword: ALL <basic>
              <HELP>
              Requests that all orbitals should be localised.
              </HELP>
              The default is
              to localise the occupied orbitals, except for PAO where the default is
              virtual.
              </KEYWORD>

  .. xmldoc:: </SELECT>

:kword:`PIPEk-Mezey`
  Requests Pipek--Mezey localisation. This is the default.

  .. xmldoc:: <SELECT MODULE="LOCALISATION" NAME="LOC_METHODS" APPEAR="Localisation method" CONTAINS="PIPE,BOYS,EDMI,CHOL,PAO,SKIP">

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="PIPE" APPEAR="Pipek-Mezey" LEVEL="ADVANCED" KIND="SINGLE" EXCLUSIVE="BOYS,EDMI,CHOL,PAO,SKIP">
              %%Keyword: PIPE <advanced>
              <HELP>
              Requests Pipek-Mezey localisation.
              </HELP>
              This is the default.
              </KEYWORD>

:kword:`BOYS`
  Requests Boys localisation. The default is Pipek--Mezey.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="BOYS" APPEAR="Boys-Forster" LEVEL="ADVANCED" KIND="SINGLE" EXCLUSIVE="PIPE,EDMI,CHOL,PAO,SKIP">
              %%Keyword: BOYS <advanced>
              <HELP>
              Requests Boys localisation.
              </HELP>
              The default is Pipek-Mezey.
              </KEYWORD>

:kword:`EDMIston-Ruedenberg`
  Requests Edmiston--Ruedenberg localisation. The default is Pipek--Mezey.
  Note that this option requires that the Cholesky (or RI/DF) representation
  of the two-electron integrals has been produced by :program:`SEWARD`.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="EDMI" APPEAR="Edmiston-Ruedenberg" LEVEL="ADVANCED" KIND="SINGLE" EXCLUSIVE="PIPE,BOYS,CHOL,PAO,SKIP">
              %%Keyword: EDMI <advanced>
              <HELP>
              Requests Edmiston-Ruedenberg localisation.
              </HELP>
              The default is Pipek-Mezey.
              Note that this option requires that the Cholesky (or RI/DF) representation
              of the two-electron integrals has been produced by SEWARD.
              </KEYWORD>

:kword:`CHOLesky`
  Requests Cholesky localisation (non-iterative). The default is Pipek--Mezey.
  This and PAO are the only options that can handle point group symmetry.
  The decomposition threshold is by default 1.0d-8 but may be changed
  through the :kword:`THREshold` keyword.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="CHOL" APPEAR="Cholesky" LEVEL="ADVANCED" KIND="SINGLE" EXCLUSIVE="PIPE,BOYS,EDMI,PAO,SKIP">
              %%Keyword: CHOL <advanced>
              <HELP>
              Requests Cholesky localisation.
              </HELP>
              The default is Pipek-Mezey.
              </KEYWORD>

:kword:`PAO`
  Requests PAO localisation (non-iterative) using Cholesky decomposition
  to remove linear dependence.
  The default is Pipek--Mezey.
  This and Cholesky are the only options that can handle point group symmetry.
  The decomposition threshold is by default 1.0d-8 but may be changed
  through the :kword:`THREshold` keyword.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="PAO" APPEAR="PAO" LEVEL="ADVANCED" KIND="SINGLE" EXCLUSIVE="PIPE,BOYS,EDMI,CHOL,SKIP">
              %%Keyword: PAO <advanced>
              <HELP>
              Requests PAO localisation.
              </HELP>
              The default is Pipek-Mezey.
              </KEYWORD>

:kword:`SKIP`
  Leaves the input orbitals unchanged. It is turned off by default.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="SKIP" APPEAR="None" LEVEL="ADVANCED" KIND="SINGLE" EXCLUSIVE="PIPE,BOYS,EDMI,CHOL,PAO">
              %%Keyword: SKIP <advanced>
              <HELP>
              Leaves the input orbitals unchanged.
              </HELP>
              </KEYWORD>

  .. xmldoc:: </SELECT>

:kword:`ITERations`
  The following line specifies the maximum number of iterations to be
  used by the iterative localisation procedures. The default is 300.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="ITER" APPEAR="Iterations" LEVEL="ADVANCED" KIND="INT">
              <HELP>
              Please, specify the maximum number of iterations to be
              used by the iterative localisation procedures. The default is 300.
              </HELP>
              %%Keyword: ITER <advanced>
              The following line specifies the maximum number of iterations to be
              used by the iterative localisation procedures. The default is 100.
              </KEYWORD>

:kword:`THREshold`
  The following line specifies the convergence threshold used for
  changes in the localisation functional. The default is 1.0d-6.
  For Cholesky and PAO methods, it is the decomposition threshold and
  the default is 1.0d-8.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="THRE" APPEAR="Functional threshold" LEVEL="ADVANCED" KIND="REAL">
              <HELP>
              Please, specify the convergence threshold used for
              changes in the localisation functional (default: 1.0d-6)
              or the decomposition threshold (default: 1.0d-8).
              </HELP>
              %%Keyword: THRE <advanced>
              The following line specifies the convergence threshold used for
              changes in the localisation functional. The default is 1.0d-6.
              For Cholesky and PAO methods, it is the decomposition threshold and
              the default is 1.0d-8.
              </KEYWORD>

:kword:`THRGradient`
  The following line specifies the convergence threshold used for
  the gradient of the localisation functional. The default is 1.0d-2.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="THRG" APPEAR="Gradient threshold" LEVEL="ADVANCED" KIND="REAL">
              <HELP>
              Please, specify the convergence threshold used for
              changes in the gradient of the localisation functional. The default is 1.0d-2.
              </HELP>
              %%Keyword: THRG <advanced>
              The following line specifies the convergence threshold used for
              the gradient of the localisation functional. The default is 1.0d-2.
              </KEYWORD>

:kword:`THRRotations`
  The following line specifies the screening threshold used in
  the Jacobi sweep optimization algorithm. The default is 1.0d-10.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="THRR" APPEAR="Screening threshold" LEVEL="ADVANCED" KIND="REAL">
              <HELP>
              Please, specify the convergence threshold used in
              the Jacobi sweep optimization algorithm. The default is 1.0d-10.
              </HELP>
              %%Keyword: THRR <advanced>
              The following line specifies the screening threshold used in
              the Jacobi sweep optimization algorithm. The default is 1.0d-10.
              </KEYWORD>

:kword:`CHOStart`
  Requests that iterative localisation procedures use Cholesky orbitals
  as initial orbitals. The default is to use the orbitals from
  :file:`INPORB` directly.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="CHOS" APPEAR="Cholesky guess" LEVEL="ADVANCED" KIND="SINGLE">
              %%Keyword: CHOS <advanced>
              <HELP>
              Requests that the localisation procedure uses Cholesky orbitals
              as initial orbitals.
              </HELP>
              The default is not to use Cholesky orbitals.
              </KEYWORD>

:kword:`ORDEr`
  Requests that the localised orbitals are ordered in the same way
  as the Cholesky orbitals would be. This is mainly useful when
  comparing orbitals from different localisation schemes. The
  ordering is done according to maximum overlap with the
  Cholesky orbitals. The default is not to order.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="ORDE" APPEAR="Orbital reordering" LEVEL="ADVANCED" KIND="SINGLE">
              %%Keyword: ORDE <advanced>
              <HELP>
              Requests that the localised orbitals are ordered in the same way
              as the Cholesky orbitals would be.
              </HELP>
              The default is not to order.
              </KEYWORD>

:kword:`DOMAin`
  Requests orbital domains and pair domains are set up and analyzed.
  The default is not to set up domains.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="DOMA" APPEAR="Orbital and pair domains analysis" LEVEL="ADVANCED" KIND="SINGLE">
              %%Keyword: DOMA <advanced>
              <HELP>
              Requests orbital domains and pair domains are set up and analyzed.
              </HELP>
              </KEYWORD>
              The default is not to set up domains.

:kword:`THRDomain`
  The following line specifies two thresholds to be used in defining
  orbital domains. The first is the Mulliken population threshold
  such that atoms are included in the domain until the population
  (divided by 2) is larger than this number (default: 9.0d-1).
  The second threshold is used for the Pulay completeness check of
  the domain (default: 2.0d-2).

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="THRD" APPEAR="Domain thresholds" LEVEL="ADVANCED" KIND="REALS" SIZE="2" REQUIRE="DOMA">
              <HELP>
              Please, specify two thresholds:
              The first is the Mulliken population threshold
              such that atoms are included in the domain (default: 9.0d-1).
              The second threshold is used for the Pulay completeness check of
              the domain (default: 2.0d-2).
              </HELP>
              </KEYWORD>
              %%Keyword: THRD <advanced>
              The following line specifies two thresholds to be used in defining
              orbital domains. The first is the Mulliken population threshold
              such that atoms are included in the domain until the population
              (divided by 2) is larger than this number (default: 9.0d-1).
              The second threshold is used for the Pulay completeness check of
              the domain (default: 2.0d-2).

:kword:`THRPairdomain`
  The following line specifies three thresholds to be used for
  classifying pair domains: R1, R2, and R3. (Defaults: 1.0d-10,
  1.0d1, and 1.5d1.)
  If R is the smallest distance
  between two atoms in the pair domain (union of the individual orbital
  domains), then pair domains are classified according to:
  R :math:`\leq` R1: strong pair,
  R1 :math:`<` R :math:`\leq` R2: weak pair,
  R2 :math:`<` R :math:`\leq` R3: distant pair, and
  R3 :math:`<` R: very distant pair.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="THRP" APPEAR="Pair domain threshold" LEVEL="ADVANCED" KIND="REALS" SIZE="3" REQUIRE="DOMA">
              <HELP>
              Please, specify three thresholds to be used for
              classifying pair domains: R1, R2, and R3. (Defaults: 1.0d-10,
              1.0d1, and 1.5d1.)
              </HELP>
              </KEYWORD>
              %%Keyword: THRP <advanced>
              The following line specifies three thresholds to be used for
              classifying pair domains: R1, R2, and R3. (Defaults: 1.0d-10,
              1.0d1, and 1.5d1.)

:kword:`LOCNatural orbitals`
  This keyword is used to select atoms for defining the localised natural
  orbitals (LNOs), thus a set of localised orbitals with well-defined occupation numbers.
  All other options specified in the :program:`LOCALISATION` program input apply (e.g., input orbitals,
  localisation method, etc.).
  On the next line give the number of atoms that identify the region of interest
  and the threshold used to select the localised orbitals belonging to this region
  (recommended values > 0.2 and < 1).
  An additional line gives the names of the (symmetry unique) atoms as defined in the :program:`SEWARD` input.
  The keyword :kword:`LOCN` is used to define suitable occupation numbers for RASSCF active orbitals
  that have been localised. It has proven useful in Effective Bond Order (EBO) analysis.
  Here is a sample input for a complex containing an iron-iron multiple bond. ::

    LOCN
    2  0.3
    Fe1  Fe2

  In this example, the (localised) orbitals constructed by the :program:`LOCALISATION` program
  are subdivided in two groups: those having less than 0.3 total Mulliken population on
  the two iron atoms, and the remaining orbitals, obviously localised on the iron-iron region. The resulting
  density matrices for the two subsets of orbitals are then diagonalized separately
  and the corresponding (localised) natural orbitals written to :file:`LOCORB` with the proper occupation
  numbers. Note that the two sets of LNOs are mutually non-orthogonal.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="LOCN" APPEAR="Localized natural orbitals" LEVEL="BASIC" KIND="CUSTOM">
              <HELP>
              Specify the number of atoms in the region and the threshold.
              Then the names of the symmetry unique atoms.
              </HELP>
              </KEYWORD>
              %%Keyword: LOCN <basic>
              This keyword is used to select atoms for defining the localised natural
              orbitals (LNOs), thus a set of localised orbitals with well-defined occupation numbers.
              All other options specified in the localisation input apply (e.g., input orbitals,
              localisation method, etc.).
              On the next line give the number of (symmetry unique) atoms that identify the region of interest
              and the threshold used to select the localised orbitals belonging to this region.
              An additional line gives the names of the atoms as defined in the SEWARD input.
              This keyword is used to define occupation numbers when localising active orbitals
              from RASSCF calculations. Particularly useful in Effective Bond Order (EBO) analysis.

:kword:`LOCCanonical orbitals`
  This keyword is used to select atoms for defining the localised canonical
  orbitals (LCOs), thus a set of localised orbitals with well-defined orbital energies
  (eigenvalues of a local Fock matrix).
  Please, refer to the analogous keyword :kword:`LOCN` in this manual for more details and input examples.

  .. xmldoc:: <KEYWORD MODULE="LOCALISATION" NAME="LOCC" APPEAR="Localized canonical orbitals" LEVEL="BASIC" KIND="CUSTOM">
              <HELP>
              Specify the number of atoms in the region and the threshold.
              Then the names of the symmetry unique atoms.
              </HELP>
              </KEYWORD>
              %%Keyword: LOCC <basic>
              This keyword is used to select atoms for defining the localised canonical
              orbitals (LCOs), thus a set of localised orbitals with well-defined orbital energies.
              All other options specified in the localisation input apply (e.g., input orbitals,
              localisation method, etc.).
              On the next line give the number of (symmetry unique) atoms that identify the region of interest
              and the threshold used to select the localised orbitals belonging to this region.
              An additional line gives the names of the atoms as defined in the SEWARD input.

Limitations
...........

The limitations on the number of basis functions are the same as specified
for :program:`SEWARD`.

Input examples
..............

This input is an example of the Boys localisation of the CO molecule. Note that no
symmetry should not be used in any calculation of localised orbitals except for
Cholesky and PAO orbitals.

.. extractfile:: ug/localisation.Boys.input

  &GATEWAY
  Coord = $MOLCAS/Coord/CO.xyz
  Basis = STO-3G
  Group = C1

  &SEWARD ; &SCF

  &LOCALISATION
  Boys

This input is an example of the Projected Atomic Orbital localisation of the
virtual orbitals of the CO molecule. The threshold for the Cholesky
decomposition that removes linear dependence is set to 1.0d-14.

.. extractfile:: ug/localisation.PAO.input

  &GATEWAY
  Coord = $MOLCAS/Coord/CO.xyz
  Basis = STO-3G
  Group = C1

  &SEWARD ; &SCF

  &LOCALISATION
  PAO
  Threshold = 1.0d-14

This input is an example of the Cholesky localisation (using default 1.0d-8 as
threshold for the decomposition) of the
valence occupied orbitals of the CO molecule.
Orbital domains are set up and analyzed.

.. extractfile:: ug/localisation.Cholesky.input

  &GATEWAY
  Coord = $MOLCAS/Coord/CO.xyz
  Basis = STO-3G
  Group = C1

  &SEWARD ; &SCF

  &LOCALISATION
  Cholesky
  Freeze
  Domain

.. xmldoc:: </MODULE>
