.. index::
   single: Program; SCF
   single: SCF

.. _UG\:sec\:scf:

:program:`scf`
==============

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. _UG\:sec\:scf_description:

Description
-----------

.. xmldoc:: <MODULE NAME="SCF">
            %%Description:
            <HELP>
            The SCF program of the molcas program system generates
            closed-shell Hartree-Fock, open-shell UHF, and Kohn Sham DFT wave functions.
            </HELP>

The :program:`SCF` program of the |molcas| program system generates
closed-shell Hartree--Fock, open-shell UHF, and Kohn Sham DFT wave functions.

The construction of the Fock
matrices is either done conventionally from the two-electron integral
file :file:`ORDINT`,
which was generated in a previous step by the :program:`SEWARD`
integral code, or alternatively (only for closed shell calculations)
integral-direct by recomputing all the
two-electron integrals when needed :cite:`AlmlofFaegriKorsell_DirSCF`.
The later route is recommended for
large basis sets or molecules, when the two-electron integral file would
become extensively large. It is automatically taken, when the :program:`SCF`
program cannot find any :file:`ORDINT` file in the work directory.
The direct Fock matrix construction employs an efficient integral prescreening
scheme, which is based on differential densities
:cite:`CremerGauss_SCFprescreening,HaeserAhlrichs_SCFprescreening`: only
those AO integrals are computed, where the estimated contractions with the
related differential density matrix elements give significant (Coulomb or
exchange) contributions to the (differential) two-electron part of the Fock
matrix. Integral prescreening is performed at two levels, (i) at the level
of shell quadruples, and (ii) at the level of individual primitive Gaussians.
Prescreening at the level of contracted functions is not supported, because
this would be inefficient in the context of a general contraction scheme.
In order to work with differential density and Fock matrices, a history of
these entities over previous iterations has to be kept. All these matrices
are partly kept in memory, and partly held on disk. The :program:`SCF` program
either works with simple differences of the actual and the previous density,
or alternatively with minimized densities, obtained from linear combinations
of the actual density and all the previous minimized densities.

Besides the conventional and the fully-direct algorithms there is also a
semi-direct path, which allows for the storage of some of the AO integrals
during the first iteration, which then are retrieved from disk in subsequent
iterations. That path is taken, if the keyword :kword:`DISK` with an
appropriate
argument specifying the amount of AO integrals to store is found on the
input stream. The semi-direct path is recommended for medium sized problems,
where the two-electron integral file would become a bit too large (but not
orders of magnitude).

.. compound::

  The program contains a feature that allows you to make the
  orbitals partially populated during the aufbau procedure.
  This feature is not primarily intended to accelerate the convergence
  but rather to ensure that you do get convergence in
  difficult cases.
  The orbitals are populated with with electrons according to

  .. math:: \eta_i=\frac{2}{1+e^{(\varepsilon_i-\varepsilon_f)/kT}}

  where
  :math:`\varepsilon_i`
  is the orbital energy of orbital :math:`i` and
  :math:`\varepsilon_f`
  is the Fermi energy.
  In this "Fermi aufbau" procedure
  the temperature is slowly lowered until it reaches a minimum
  value and then kept constant until
  a stable closed shell configuration is determined.
  Then normal SCF iterations will be performed with the selected
  closed shell configuration.
  For systems that are not really closed shell systems, for example
  diradicals, you might end up in the situation that the program
  does not find any stable closed shell configuration.
  In that case it will continue to optimize the closed shell
  energy functional with partial occupation numbers.
  If this is the case, this is probably what you want, and such
  orbitals would be ideal as starting orbitals for an MCSCF
  calculation.

The initial orbital guess is either obtained by diagonalizing the bare nuclei
Hamiltonian, from an initial guess produced by the module :program:`Guessorb`
or from orbitals of a previous Hartree--Fock SCF calculation.
These starting orbitals are automatically located in the order

#. SCF orbitals from a previous calculation located in the :file:`RUNFILE`

#. SCF orbitals from a previous calculation located in a formatted orbitals file, :file:`INPORB`.

#. initial guess orbitals from module :program:`Guessorb` located in the :file:`RUNFILE` and

The program has three types of convergence accelerating schemes:
(i) dynamic damping :cite:`Karlstroem_SCFdyndumping`, (ii) the :math:`C^2`\-DIIS method
using the orbital gradient as error vector :cite:`c2-diis`, and (iii)
a combined second-order update/\ :math:`C^2`\-DIIS procedure. The latter eliminates the
Brillouin violating elements of the Fock matrix by proper orbital rotations
and hence avoids diagonalization of the Fock matrix: the approximate inverse
Hessian is updated (BFGS) in a first step, and then the new orbital
displacement vector is obtained from the updated Hessian using :math:`C^2`\-DIIS
extrapolation :cite:`FischerAlmloef_OrbRot`.
Dynamic damping gives substantial improvements in highly anharmonic
regions of the energy hyper surface, while the second-order
update/\ :math:`C^2`\-DIIS procedure exhibits excellent convergence for less anharmonic
regions. By default, dynamic damping is used during the first few iterations.
When the change in the density between two subsequent iterations drops below
a certain threshold the second-order update/\ :math:`C^2`\-DIIS procedure kicks in.
It is also possible to use the older first order :math:`C^2`\-DIIS scheme instead of
the second-order update/\ :math:`C^2`\-DIIS procedure by setting the density threshold
for the latter to zero in the corresponding input card (keyword
:kword:`QNRThreshold`).

By default :program:`SCF` behaves in different ways depending on what
kind of start orbitals are found according to

#. No start orbitals are found. In this case the core hamiltonian
   is diagonalized and these orbitals are used as start.
   The "Fermi aufbau" procedure is used until a stable configuration is found.

#. Start orbitals from :program:`Guessorb` are found.
   In this case the HOMO LUMO gap is analyzed and if it is small
   the "Fermi aufbau" procedure is used until a stable configuration is found.
   Otherwise the configuration suggested by :program:`Guessorb` is used.

#. Start orbitals from a previous :program:`SCF` calculation is found.
   The configuration from the previous :program:`SCF` calculation is used,
   unless some problem is detected such as partial occupation numbers
   from an unconverged calculation. In the latter case "Fermi aufbau" is used.

#. Start orbitals from an :file:`INPORB` is in the same way as for
   start orbitals from an :program:`SCF` calculation, see last point.

This behavior can be changed by suitable keywords described below.

One of the main objects of the :program:`SCF` program in the context of the
|molcas| program system is to generate starting
orbitals for subsequent MCSCF calculations.
Two options are available to
improve the canonical Hartree--Fock orbitals in this respect.

(i) It is possible to specify pseudo occupation numbers that are neither
zero nor two, thus simulating to some extent an open shell system. The
resulting wavefunction does not have any physical meaning, but will
provide better starting orbitals for open shell systems.

(ii) Usually, the lowest virtual canonical Hartree--Fock orbitals are
too diffuse as correlating orbitals in an MCSCF calculation.
If the keyword :kword:`IVO` is encountered in the input stream, the
:program:`SCF` program will diagonalize the core Hamiltonian matrix within
the virtual space and write the resulting more compact eigenvectors to the
:file:`SCFORB` and :file:`RUNFILE` files,
rather than the virtual eigenvectors of the Fock
matrix. It should be noted, that this option must never be used, if the
:program:`SCF` wave function itself is used subsequently as a reference
function: no MP2 or coupled cluster calculations after an
:program:`SCF` run with :kword:`IVO`!

A further method to generate starting orbitals for MCSCF calculations is
to perform an SCF calculation for a slightly positively charged moiety.

.. _UG\:sec\:scf_dependencies:

Dependencies
------------

The :program:`SCF` program requires the one-electron integral file
:file:`ONEINT` and the communications file :file:`RUNFILE`,
which contains among others the
basis set specifications processed by :program:`SEWARD`. For conventional
(not integral-direct) runs the two-electron integral file :file:`ORDINT`
is required as well. All these files are generated by a preceding
:program:`SEWARD` run.

.. index::
   pair: Files; SCF

.. _UG\:sec\:scf_files:

Files
-----

Below is a list of the files that are used/created by the program
:program:`SCF`.

Input files
...........

:program:`SCF` will use the following input
files: :file:`ONEINT`, :file:`ORDINT`, :file:`RUNFILE`, :file:`INPORB`
(for more information see :ref:`UG:sec:files_list`).

Output files
............

.. class:: filelist:

:file:`SCFORB`
  :program:`SCF` orbital output file.
  Contains the canonical Hartree--Fock orbitals for closed shell calculations.
  If the :kword:`IVO` option
  was specified, the virtual orbitals instead are those that diagonalize the bare
  nuclei Hamiltonian within that subspace.

:file:`UHFORB`
  Contains the canonical Hartree--Fock orbitals for open shell calculations.

:file:`UNAORB`
  This file is produced if you make a UHF calculation and it contain
  natural orbitals.

:file:`MD_SCF`
  Molden input file for molecular orbital analysis.

.. index::
   pair: Input; SCF

.. _UG\:sec\:scf_inpscf:

Input
-----

Below follows a description of the input to :program:`SCF`.

The input for each module is preceded by its name like: ::

  &SCF

Argument(s) to a keyword, either individual or composed by several entries,
can be placed in a separated line or in the same line separated by a semicolon.
If in the same line, the first argument requires an equal sign after the
name of the keyword.

Basic general keywords
......................

Below is a list of keywords that should cover the needs of most users.

.. class:: keywordlist

:kword:`TITLe`
  One line for the title

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="TITLE" KIND="STRING" LEVEL="BASIC">
              %%Keyword: Title <basic>
              <HELP>
              One line for the title
              </HELP>
              </KEYWORD>

:kword:`UHF`
  Use this keyword to run Unrestricted Hartree--Fock code.
  Note that current implementation of UHF code has some
  restrictions, and not all features of :program:`SCF` program are supported.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="UHF" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: UHF <basic> GUI:keyword
              <HELP>
              Use this keyword to run Unrestricted Hartree-Fock code
              Note that current implementation of UHF code has some
              restrictions, and not all features of SCF program are supported
              </HELP>
              </KEYWORD>

:kword:`ZSPIN`
  Use this keyword to specify the difference in the number of :math:`\alpha` and :math:`\beta`
  electrons in the system. The default is 0 or 1 depending on if there is an even
  or odd number of electrons.
  Any value different from 0 requires the :kword:`UHF` keyword.
  This keyword is not needed when you specify the number of electrons with
  the keyword :kword:`OCCUpied`.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="ZSPIN" APPEAR="zSpin" KIND="INT" LEVEL="BASIC" REQUIRE="UHF" EXCLUSIVE="SPIN">
              %%Keyword: zSpin <basic>
              <HELP>
              Use this keyword to specify the difference in the number of alpha
              and beta electrons in the system. The default is 0 or 1 depending
              on if there is an even or odd number of electrons. This keyword
              is not needed when you specify the number of electrons with the
              keyword OCCUpied.
              </HELP>
              </KEYWORD>

:kword:`SPIN`
  Alternative way of specifying the electronic spin of the system.
  The keyword is followed by an integer giving the value of spin multiplicity (:math:`2S+1`).
  Default is 1 (singlet) or 2 (doublet) depending on if there is an even or odd number of electrons.
  Any value different from 1 requires the :kword:`UHF` keyword.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="SPIN" APPEAR="Spin" LEVEL="BASIC" KIND="INT" MIN_VALUE="1" REQUIRE="UHF" EXCLUSIVE="ZSPIN">
              %%Keyword: Spin <basic>
              <HELP>
              The keyword is followed by an integer giving the value of spin
              multiplicity (2S+1). Default is 1 (singlet) or 2 (doublet)
              depending on if there is an even or odd number of electrons.
              </HELP>
              </KEYWORD>

:kword:`KSDFT`
  Use this Keyword to do the Density Functional Theory Calculations.
  This Keyword should be followed by functional Keyword:
  BLYP, B3LYP, B3LYP5, HFB, HFS, LDA, LDA5, LSDA, LSDA5, SVWN, SVWN5, TLYP, PBE, PBE0, M06, M06HF, M062X, M06L.
  Example: KSDFT=B3LYP

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="KSDFT" APPEAR="DFT" KIND="CHOICE" LIST="----,BLYP,B3LYP,B3LYP5,HFB,HFS,LDA,LDA5,LSDA,LSDA5,SVWN,SVWN5,TLYP,PBE,PBE0,M06,M06HF,M062X,M06L" LEVEL="BASIC">
              %%Keyword: KSDFT <basic>
              <HELP>
              Use this Keyword to do the Density Functional Theory Calculations
              This Keyword should be followed by the functional Keyword:
              BLYP, B3LYP, B3LYP5, HFB, HFS, LDA, LDA5, LSDA, LSDA5, SVWN, SVWN5, TLYP, PBE, PBE0, M06, M062X, M06HF, M06L.
              Example: KSDFT=B3LYP
              </HELP>
              </KEYWORD>

:kword:`DFCF`
  Use this Keyword to scaled the exchange terms and/or correlation terms of a density functional.
  This Keyword should be followed by the scaling factor for the exchange terms and the scaling factor for the correlation terms, separated by a space.
  If the values are 1.0 (default), then the original density functional is used. 
  For an HLE-type functional, use 1.25 (for exchange) and 0.5 (for correlation).
  Example: DFCF=1.25 0.5

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="DFCF" APPEAR="DFT exch. & corr. scaling factors" KIND="REALS" SIZE="2" LEVEL="ADVANCED">
              %%Keyword: DFCF <advanced>
              <HELP>
              Use this Keyword to scaled the exchange terms and/or correlation terms of a density functional.
              This Keyword should be followed by the scaling factor for the exchange terms 
              and the scaling factor for the correlation terms, separated by a space.
              If the values are 1.0 (default), then the original density functional is used.
              For an HLE-type functional, use 1.25 (for exchange) and 0.5 (for correlation).
              Example: DFCF=1.25 0.5
              </HELP>
              </KEYWORD>

:kword:`CHARge`
  .. compound::

    Use this keyword to set the number of electrons in the system.
    This number is defined by giving the net charge of the system.
    If this keyword is not specified, the molecule is assumed to
    have net charge zero.
    The input is given as ::

      Charge=n

    where ``n`` is the charge of the system.

  .. xmldoc:: <SELECT MODULE="SCF" NAME="ELECTRONS" APPEAR="Electron Count" LEVEL="BASIC" CONTAINS="DEFAULT,CHARGE,OCCUPIED">

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="CHARGE" KIND="INT" LEVEL="BASIC" EXCLUSIVE="OCCUPIED" MEMBER="ELECTRONS">
              %%Keyword: Charge <basic>
              <HELP>
              Use this keyword to set the number of electrons in the system.
              This number is defined by giving the net charge of the system.
              If this keyword is not specified, the molecule is assumed to
              have net charge zero.
              The input is given as
              ||
              ||Charge=n
              ||
              where n is the charge of the system.
              </HELP>
              </KEYWORD>

:kword:`OCCUpied`
  .. compound::

    Use this keyword to set the number of electrons in the system.
    This number is defined by giving the number of electron pairs
    per irreducible representation of the subgroup of :math:`D_{2h}` used
    in the calculation.
    You can use one and only one of the keywords,
    :kword:`CHARge` and
    :kword:`OCCUpied` for this purpose.
    If neither of these keywords are specified
    :kword:`CHARge` is assumed with a net charge of zero.
    It should be noted that the "Fermi aufbau"
    procedure is not used when you specify this keyword.
    The input for one of the point groups :math:`D_2`, :math:`C_{2h}` or :math:`C_{2v}`
    is given as ::

      OCCUpied= n1 n2 n3 n4

    where ``n1`` is the number of electron pairs (occupied orbitals)
    in the first irreducible representation, etc.

  If :kword:`UHF` keyword was specified, occupation numbers
  must be specified in two lines: for alpha and beta spins

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="OCCUPIED" KIND="INTS_LOOKUP" SIZE="NSYM" LEVEL="BASIC" EXCLUSIVE="CHARGE" MEMBER="ELECTRONS">
              %%Keyword: Occupied <basic>
              <HELP>
              Use this keyword to set the number of electrons in the system.
              This number is defined by giving the number of electron pairs
              per irreducible representation of the subgroup of D2h used in
              the calculation. You can use one and only one of the keywords,
              CHARge and OCCUpied for this purpose. If neither of
              these keywords are specified CHARge is assumed with a net
              charge of zero. It should be noted that the "fermi aufbau"
              procedure is not used when you specify this keyword.
              The input for one of the point groups D2, C2h or C2v
              is given as
              ||
              ||OCCUpied= n1 n2 n3 n4
              ||
              where n1 is the number of electron pairs (occupied orbitals)
              in the first irreducible representation, etc.
              If UHF keyword was specified, occupation numbers
              must be specified in two lines: for alpha and beta spins
              </HELP>
              </KEYWORD>

  .. xmldoc:: </SELECT>

:kword:`FERMi`
  .. compound::

    Use this keyword to specify that you want to use the "Fermi aufbau"
    procedure for the first few iterations to ensure convergence.
    The orbitals will be partially populated according to a Fermi
    population.
    The input is gives as ::

      Fermi= m

    where ``m`` is the temperature parameter according to

  * ``m=0``: No temperature is used. Not recommended.
  * ``m=1``: A low temperature is used and will yield swift convergence
    for well behaved systems.
  * ``m=2``: A medium low temperature is used and will yield swift and
    safe convergence for most systems. This is the default value.
  * ``m=3``: A medium temperature is used and you will obtain good
    convergence for closed shell systems. If the system is not
    a closed shell system, the temperature dependent aufbau
    procedure may not terminate. This will result in a density
    matrix with fractional occupation numbers.
  * ``m=4``: A medium high temperature is used and the temperature
    dependent aufbau procedure will most probably not terminate.
    This is useful for generating starting orbitals for an MCSCF
    calculation.
  * ``m=5``: A high temperature is used. Behaves as m=4 only more so.

  It should be noted that only dynamic damping is used until the
  program have found a stable closed shell configuration. When
  this have happened the more efficient methods: the ordinary
  :math:`C^2`\-DIIS and the second order update/\ :math:`C^2`\-DIIS procedure, are
  enabled.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="FERMI" KIND="INT" LEVEL="BASIC">
              %%Keyword: Fermi <basic> GUI:number
              <HELP>
              Use this keyword to specify that you want to use the "Fermi aufbau"
              procedure for the first few iterations to ensure convergence.
              The orbitals will be partially populated according to a Fermi
              population.
              The input is gives as
              ||
              ||Fermi= m
              ||
              where m is the temperature parameter according to
              ||
              ||m=0: No temperature is used. Not recommended.
              ||m=1: A low temperature is used and will yield swift convergence
              ||     for well behaved systems.
              ||m=2: A medium low temperature is used and will yield swift and
              ||     safe convergence for most systems. This is the default value.
              ||m=3: A medium temperature is used and you will obtain good
              ||     convergence for closed shell systems. If the system is not
              ||     a closed shell system, the temperature dependent aufbau
              ||     procedure may not terminate. This will result in a density
              ||     matrix with fractional occupation numbers.
              ||m=4: A medium high temperature is used and the temperature
              ||     dependent aufbau procedure will most probably not terminate.
              ||     This is useful for generating starting orbitals for an MCSCF
              ||     calculation.
              ||m=5: A high temperature is used. Behaves as m=4 only more so.
              ||
              It should be noted that only dynamic damping is used until the
              program have found a stable closed shell configuration. When
              this have happened the more efficient methods: the ordinary
              C2-DIIS and the second order update/C2-DIIS procedure, are
              enabled.
              </HELP>
              </KEYWORD>

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="ITER" APPEAR="Max Iterations" KIND="INT" LEVEL="BASIC">
              %%Keyword: Iterations <basic>
              <HELP>
              Specifies the maximum number of iterations. The default is 400 which
              is also the largest number you can specify.
              </HELP>
              </KEYWORD>

:kword:`CHOLesky`
  :program:`SCF` will use Cholesky (or RI/DF) representation of the two-electron integrals to compute
  the corresponding contributions to the Fock or KS matrices. The default (LK)
  algorithm is used. The configuration may be tailored using the ChoInput section.
  Default is to not use Cholesky unless the Cholesky (or RI/DF) representation of the two-electron
  integrals has been produced by :program:`SEWARD`.

  .. xmldoc:: %%Keyword: Cholesky <advanced>
              Use of Cholesky (or RI/DF) representation for the two-electron integrals
              with default SCF settings.

:kword:`CHOInput`
  This marks the start of an input section for modifying
  the default settings of the Cholesky SCF.
  Below follows a description of the associated options.
  The options may be given in any order,
  and they are all optional except for
  :kword:`ENDChoinput` which marks the end of the :kword:`CHOInput` section.

  * :kword:`NoLK`
    Available only within :kword:`ChoInput`. Deactivates the "Local Exchange" (LK) screening algorithm :cite:`Aquilante:07a` in computing
    the Fock matrix. The loss of speed compared to the default algorithm can be substantial, especially for electron-rich systems.
    Default is to use LK.

    .. xmldoc:: <GROUP MODULE="SCF" NAME="CHOINPUT" APPEAR="Cholesky input section" KIND="BLOCK" LEVEL="ADVANCED">
                %%Keyword: Choinput <advanced>
                <HELP>
                Manually modify the settings of the Cholesky SCF.
                </HELP>

    .. xmldoc:: <KEYWORD MODULE="SCF" NAME="NOLK" APPEAR="Turn off LK screening" LEVEL="ADVANCED" KIND="SINGLE">
                %%Keyword: NoLK <advanced>
                <HELP>
                Deactivates LK screening. Available only within ChoInput.
                </HELP>
                </KEYWORD>

  * :kword:`DMPK`
    Available only within :kword:`ChoInput`. Modifies the thresholds used in the LK screening.
    The keyword takes as argument a (double precision) floating point (non-negative) number used
    as correction factor for the LK screening thresholds.
    The default value is 1.0d0. A smaller value results in a slower but more accurate calculation.

    **Note:** The default choice of the LK screening thresholds is tailored to achieve as much as possible an
    accuracy of the converged SCF energy consistent with the choice of the Cholesky decomposition
    threshold.

    .. xmldoc:: <KEYWORD MODULE="SCF" NAME="DMPK" APPEAR="Damping for LK" LEVEL="ADVANCED" KIND="REAL" EXCLUSIVE="NOLK">
                %%Keyword: dmpK <advanced>
                <HELP>
                Modifies the thresholds used in the LK screening. Available only within ChoInput
                The default value is 1.0d0. A smaller value results in a slower but more accurate calculation.
                </HELP>
                </KEYWORD>

  * :kword:`NODEcomposition`
    Available only within :kword:`ChoInput`. Deactivates the Cholesky decomposition of the AO 1-particle density matrix.
    The Exchange contribution to the Fock matrix is therefore computed using occupied canonical orbitals
    instead of (localized) "Cholesky MOs" :cite:`Aquilante:06a`. This choice tends to lower the performances of the
    LK screening.
    Default is to perform this decomposition in order to obtain the Cholesky MOs.

    .. xmldoc:: <KEYWORD MODULE="SCF" NAME="NODE" APPEAR="Turn off density decomposition" LEVEL="ADVANCED" KIND="SINGLE">
                %%Keyword: NODE <advanced>
                <HELP>
                The Exchange contribution to the Fock matrix is computed using occupied canonical orbitals
                instead of (localized) "Cholesky MOs". Available only within ChoInput.
                </HELP>
                </KEYWORD>

    .. xmldoc:: </GROUP>

  * :kword:`TIME`
    Activates printing of the timings of each task of the Fock matrix build.
    Default is to not show these timings.

  * :kword:`MEMFraction`
    Set the fraction of memory to use as global Cholesky vector buffer.
    Default: for serial runs 0.0d0; for parallel runs 0.3d0.

:kword:`CONStraints`
  Performs a Constrained (Natural Orbitals) SCF calculation, available only in combination with Cholesky or RI integral representation.
  An example of input for the keyword :kword:`CONS` is the following: ::

    CONStraints
     2  3
     1 -1
     1  1  1

    ADDCorrelation
    pbe

    SAVErage

  The keyword :kword:`CONS` has two compulsory arguments: the number of constrained NOs
  (in each irrep) to be used in the CNO-SCF calculation, followed by one line per irrep specifying the spin configuration
  of the so-called (+) wavelet (-1 :math:`\rightarrow` beta, 1 :math:`\rightarrow` alpha)
  The OPTIONAL keyword :kword:`ADDC` is used to include a correlation energy correction through a DFT functional specified as argument (LDA, LDA5, PBE and BLYP available at the moment)
  The OPTIONAL keyword :kword:`SAVE` forces the program to use spin-averaged wavelets.

  .. xmldoc:: %%Keyword: CONS <advanced>
              Performs a Constrained (Natural Orbitals) SCF calculation, available only in combination with Cholesky or RI integral representation.
              An example of input for the keyword CONS is the following:
              ||
              ||CONStraints
              || 2  3
              || 1 -1
              || 1  1  1
              ||
              ||ADDCorrelation
              ||pbe
              ||
              ||SAVErage
              ||
              The keyword CONS has two compulsory arguments: the number of constrained NOs
              (in each irrep) to be used in the CNO-SCF calculation, followed by one line per irrep specifying the spin configuration
              of the so-called (+) wavelet (-1 --> beta, 1 --> alpha)
              The OPTIONAL keyword ADDC is used to include a correlation energy correction through a DFT functional specified as argument (LDA, LDA5,
              PBE and BLYP available at the moment)
              The OPTIONAL keyword SAVE forces the program to use spin-averaged wavelets.

:kword:`OFEMbedding`
  Performs a Orbital-Free Embedding (OFE)SCF calculation, available only in combination with Cholesky or RI integral representation.
  The runfile of the environment subsystem renamed AUXRFIL is required.
  An example of input for the keyword :kword:`OFEM` is the following: ::

    OFEMbedding
     ldtf/pbe
    dFMD
     1.0   1.0d2
    FTHAw
     1.0d-4

  The keyword :kword:`OFEM` requires the specification of two functionals in the form fun1/fun2, where fun1 is the functional
  used for the Kinetic Energy (available functionals: Thomas--Fermi, with acronym LDTF, and the NDSD functional), and where
  fun2 is the xc-functional (LDA, LDA5, PBE and BLYP available at the moment).
  The OPTIONAL keyword :kword:`dFMD` has two arguments: first, the fraction of correlation potential to be added to the
  OFE potential (zero for KSDFT and one for HF); second, the exponential decay factor for this correction (used in PES calculations).
  The OPTIONAL keyword :kword:`FTHA` is used in a freeze-and-thaw cycle (EMIL Do While) to specify the (subsystems) energy convergence threshold.

  .. xmldoc:: %%Keyword: OFEM <advanced>
              Performs a Orbital-Free Embedding (OFE)SCF calculation, available only in combination with Cholesky or RI integral representation.
              The runfile of the environment subsystem renamed AUXRFIL is required.
              An example of input for the keyword OFEM is the following:
              ||
              ||OFEMbedding
              || ldtf/pbe
              ||dFMD
              || 1.0  1.0d2
              ||FTHAw
              || 1.0d-4
              ||
              The keyword OFEM requires the specification of two functionals in the form fun1/fun2, where fun1 is the functional used for the
              Kinetic Energy (available functionals: Thomas-Fermi, with acronym LDTF, and the NDSD functional), and where
              fun2 is the xc-functional (LDA, LDA5, PBE and BLYP available at the moment).
              The OPTIONAL keyword dFMD has two arguments: first, the fraction of correlation potential to be added to the OFE potential (zero for
              KSDFT and one for HF); second, the exponential decay factor for this correction (used in PES calculations).
              The OPTIONAL keyword FTHA is used in a freeze-and-thaw cycle (EMIL Do While) to specify the (subsystems) energy convergence
              threshold.

:kword:`ITERations`
  Specifies the maximum number of iterations. The default is 400 which
  is also the largest number you can specify.

:kword:`CORE`
  The starting vectors are obtained from a diagonalization of the core
  Hamiltonian.

  .. xmldoc:: <SELECT MODULE="SCF" NAME="GUESS" APPEAR="Initial guess" LEVEL="BASIC" CONTAINS="CORE,LUMORB,FILEORB,GSSRUNFILE">

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="CORE" APPEAR="Core" KIND="SINGLE" LEVEL="BASIC" EXCLUSIVE="LUMORB,FILEORB,GSSRUNFILE" MEMBER="GUESS">
              %%Keyword: Core <basic>
              <HELP>
              The starting vectors are obtained from a diagonalization of the core
              Hamiltonian.
              </HELP>
              </KEYWORD>

:kword:`LUMORB`
  The starting vectors are taken from a previous :file:`SCFORB` file called
  :file:`INPORB`.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="LUMORB" APPEAR="LUMORB" KIND="SINGLE" LEVEL="BASIC" EXCLUSIVE="FILEORB,CORE,GSSRUNFILE" MEMBER="GUESS">
              %%Keyword: LUMORB <basic>
              <HELP>
              The starting vectors are taken from a previous SCFORB file called
              INPORB.
              </HELP>
              </KEYWORD>

:kword:`FILEORB`
  The starting vectors are taken from a previous :file:`SCFORB` file, specified by user.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="FILEORB" APPEAR="FILEORB" KIND="STRING" LEVEL="BASIC" EXCLUSIVE="LUMORB,CORE,GSSRUNFILE" MEMBER="GUESS">
              %%Keyword: FILEORB <basic>
              <HELP>
              The starting vectors are taken from a previous SCFORB file, specified by user.
              </HELP>
              </KEYWORD>

:kword:`GSSRunfile`
  The starting vectors are taken from the orbitals produced by :program:`Guessorb`.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="GSSRUNFILE" APPEAR="GSSRUNFILE" KIND="SINGLE" LEVEL="BASIC" EXCLUSIVE="CORE,LUMORB,FILEORB" MEMBER="GUESS">
              %%Keyword: GssRunfile <basic>
              <HELP>
              The starting vectors are taken from the orbitals produced by Guessorb.
              </HELP>
              </KEYWORD>

  .. xmldoc:: </SELECT>

:kword:`HLGAp`
  This keyword is used to make the program level shift the virtual
  orbitals in such a way that the HOMO LUMO gap is at least the value
  specified on the next line. This will help convergence in difficult
  cases but may lead to that it converges to an excited configuration.
  A suitable value is 0.2.

  .. xmldoc:: %%Keyword: HLgap <basic>
              This keyword is used to make the program levelshift the virtual
              orbitals in such a way that the HOMO LUMO gap is at least the value
              specified on the next line. This will help convergence in difficult
              cases but may lead to that it converges to an excited configuration.
              A suitable value is 0.2.

Advanced general keywords
.........................

.. class:: keywordlist

:kword:`SCRAmble`
  This keyword will make the start orbitals slightly scrambled,
  accomplished by making a few small random orbital rotations.
  How much the orbitals are scrambled is determined by the
  parameter read on the next entry. A reasonable choice for
  this parameter is 0.2 which correspond to maximum rotation angle
  of :math:`\arcsin 0.2`.
  Using this keyword may be useful for UHF calculations with
  same number of :math:`\alpha` and :math:`\beta` electrons that are not
  closed shell cases.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="SCRAMBLE" KIND="REAL" LEVEL="ADVANCED" REQUIRE="UHF">
              %%Keyword: Scramble <advanced>
              <HELP>
              This keyword will make the start orbitals slightly scrambled,
              accomplished by making a few small random orbital rotations.
              How much the orbitals are scrambled is determined by the
              parameter read on the next entry. A reasonable choice for
              this parameter is 0.2 which correspond to maximum rotation angle
              of arcsin(0.2).
              Using this keyword may be useful for UHF calculations with
              same number of alpha and beta electrons that are not
              closed shell cases.
              </HELP>
              </KEYWORD>

:kword:`ORBItals`
  Specifies the number of orbitals in the subspace of the full
  orbital space defined by the basis set, in which the SCF energy
  functional is optimized. The size of this subspace is given for each
  of the irreducible representations of the subgroup of :math:`D_{2h}`.
  If this keyword is not specified when starting orbitals are read, the
  full orbital space is assumed.
  The keyword takes as argument *nIrrep* (# of irreps) integers.
  **Note** that this keyword is only meaningful when the :program:`SCF`
  program is fed with input orbitals (cf. :kword:`LUMORB`).

  .. xmldoc:: %%Keyword: Orbitals <advanced>
              Specifies the number of orbitals in the subspace of the full
              orbital space defined by the basis set, in which the SCF energy
              functional is optimized. The size of this subspace is given for each
              of the irreducible representations of the subgroup of D2h.
              If this keyword is not specified when starting orbitals are read, the
              full orbital space is assumed.
              The keyword takes as argument nIrrep (# of irreps) integers.
              Note that this keyword is only meaningful when the SCF
              program is fed with input orbitals (cf. LUMORB).

:kword:`FROZen`
  Specifies the number of orbitals not optimized during
  iterative procedure. The size of this subspace is given for each
  of the irreducible representations of the subgroup of :math:`D_{2h}`.
  If this keyword is not specified the number of frozen orbitals
  is set to zero for each irreducible representation.
  If the starting vectors are obtained from a diagonalization of the bare
  nuclei Hamiltonian the atomic orbitals with the lowest one-electron energy
  are frozen. If molecular orbitals are read from :file:`INPORB` the frozen
  orbitals are those that are read in first in each symmetry.
  The keyword takes as argument *nIrrep* (# of irreps) integers.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="FROZEN" KIND="INTS_LOOKUP" SIZE="NSYM" LEVEL="ADVANCED">
              %%Keyword: Frozen <advanced>
              <HELP>
              Specifies the number of orbitals not optimized during
              iterative procedure. The size of this subspace is given for each
              of the irreducible representations of the subgroup of D2h.
              If this keyword is not specified the number of frozen orbitals
              is set to zero for each irreducible representation.
              If the starting vectors are obtained from a diagonalization of the bare
              nuclei Hamiltonian the atomic orbitals with the lowest one-electron energy
              are frozen. If molecular orbitals are read from the file INPORB the frozen
              orbitals are those that are read in first in each symmetry.
              The keyword takes as argument nIrrep (# of irreps) integers.
              </HELP>
              </KEYWORD>

:kword:`OVLDelete`
  Specifies the threshold for deleting near linear dependence in the
  basis set. The eigenvectors of the overlap matrix with eigenvalues
  less than that threshold are removed from the orbital subspace, and
  do not participate in the optimization procedure. The default value
  is 1.0d-5.
  The keyword takes as argument a (double precision) floating point number.
  Note that the :file:`SCFORB` file will contain the deleted orbitals as a
  complementary set to the actual SCF orbitals! In future use of this orbital
  file the complementary set should always be deleted from use.

  .. xmldoc:: %%Keyword: Ovldelete <advanced>
              Specifies the threshold for deleting near linear dependence in the
              basis set. The eigenvectors of the overlap matrix with eigenvalues
              less than that threshold are removed from the orbital subspace, and
              do not participate in the optimization procedure. The default value
              is 1.0d-5.
              The keyword takes as argument a (double precision) floating point number.
              Note that the SCFORB file will contain the deleted orbitals as a
              complemental set to the actual SCF orbitals! In future use of this orbital
              file the complemental set should always be deleted from use.

:kword:`PRORbitals`
  Specifies which orbitals are to be printed in the log file (standard output).
  The keyword takes as argument two integers.
  The possible values are:

  .. container:: list

    0 --- No orbitals printed.

    1 --- orbitals with orbital energies smaller than
    :math:`2E_{\text{HOMO}}-E_{\text{LUMO}}` are printed.

    2 --- followed by real number (ThrEne) --- orbitals with orbital
    energies smaller than ThrEne are printed.

  Default value is 1.

  Second (optional) argument specifies a format:

  .. container:: list

    0 --- No orbitals printed

    1 --- Print only one-electron energies and occupation numbers

    2 --- Short print format

    3 --- Extended print format

  Default value is 3 for small numbers of MOs and 2 for number of MOs > 256.

  .. xmldoc:: %%Keyword: Prorbitals <advanced>
              Specifies which orbitals are to be printed in the logfile (standard output).
              The keyword takes as argument two integers.
              The possible values of first argument are:
              ||0 --- No orbitals printed;
              ||1 --- orbitals with orbital energies smaller than
              ||      2E(homo)-E(lumo) are printed; and
              ||2 --- followed by real number (ThrEne) --- orbitals with orbital
              ||      energies smaller than ThrEne are printed.
              Default value is 1.
              Second (optional) argument specifies a format:
              ||0 --- No orbitals printed
              ||1 --- Print only one-electron energies and occupation numbers
              ||2 --- Short print format
              ||3 --- Extended print format

:kword:`PRLScf`
  Specifies the general print level of the calculation. An integer
  has to be supplied as argument.
  The default value, 1, is recommended for production calculations.

  .. xmldoc:: %%Keyword: prlscf <advanced>
              Specifies the general print level of the calculation. An integer
              has to be supplied as argument.
              The default value, 1, is recommended for production calculations.

  .. :kword:`PRLInt`
       .. compound ::

         Specifies the print level in individual subroutines, primarily of those
         related to the direct construction of the Fock matrix.
         The first argument, an integer :math:`n`, specifies the number of
         subroutines. Follow :math:`n` pairs of subroutine specifiers and
         corresponding print levels are given. E.g. the input segment ::

           PRLInt= 2; 3 7 9 8

         rises the print level in subroutine 3 from 5 (default) to 7, and in
         subroutine 9 from 5 to 8. This provides more extensive information
         about how the wave function is converging (subroutine 3), and statistics
         about the efficiency of the integral prescreening (subroutine 9). This
         option is certainly not used in production calculations.

:kword:`ROBU`
  Robust LDF integral representation (non-hybrid KS-DFT only).
  Requires Local Density Fitting (LDF) in SEWARD. This is the default for LDF.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="ROBU" APPEAR="Robust LDF integral representation" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: ROBU <advanced>
              <HELP>
              Robust LDF integral representation (non-hybrid KS-DFT only). Requires Local Density Fitting (LDF) in SEWARD. This is the default for LDF.
              </HELP>
              </KEYWORD>

:kword:`NR-2`
  Nonrobust LDF integral representation with 2-index integrals only (non-hybrid KS-DFT only).
  Requires Local Density Fitting (LDF) in SEWARD. Default is robust integral representation.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="NR-2" APPEAR="Nonrobust LDF integral representation with 2-index integrals only" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: NR-2 <advanced>
              <HELP>
              Nonrobust LDF integral representation with 2-index integrals only (non-hybrid KS-DFT only). Requires Local Density Fitting (LDF) in SEWARD. Default is robust integral representation.
              </HELP>
              </KEYWORD>

:kword:`NR-3`
  Nonrobust LDF integral representation with 3-index integrals only (non-hybrid KS-DFT only).
  Requires Local Density Fitting (LDF) in SEWARD. Default is robust integral representation.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="NR-3" APPEAR="Nonrobust LDF integral representation with 3-index integrals only" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: NR-3 <advanced>
              <HELP>
              Nonrobust LDF integral representation with 3-index integrals only (non-hybrid KS-DFT only). Requires Local Density Fitting (LDF) in SEWARD. Default is robust integral representation.
              </HELP>
              </KEYWORD>

:kword:`XIDI`
  Use exact integral diagonal blocks with LDF.
  Reduces the risk of negative eigenvalues of the approximate integral matrix.
  Default is to not use exact integral diagonal blocks.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="XIDI" APPEAR="Use exact integral diagonal blocks with LDF" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: XIDI <advanced>
              <HELP>
              Use exact integral diagonal blocks with LDF. Reduces the risk of negative eigenvalues of the approximate integral matrix. Default is to not use exact integral diagonal blocks.
              </HELP>
              </KEYWORD>

:kword:`THREsholds`
  Specifies convergence thresholds. Four individual thresholds are specified
  as arguments, which have to be fulfilled simultaneously to reach convergence:
  EThr, DThr and FThr
  specify the maximum permissible difference in energy, density matrix elements
  and Fock matrix elements, respectively, in the last two iterations. The
  DltNTh finally specifies the norm of the orbital displacement vector used
  for the orbital rotations in the second-order/\ :math:`C^2`\-DIIS procedure.
  The corresponding values are read in the order given above.
  The default values are 1.0d-9, 1.0d-4, 1.5d-4, and 0.2d-4,
  respectively.
  **Note** that these thresholds automatically define the threshold
  used in the direct Fock matrix construction to estimate individual
  contributions to the Fock matrix such that
  the computed energy will have an accuracy that is better than the
  convergence threshold.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="THRESHOLD" KIND="REALS" SIZE="4" LEVEL="ADVANCED">
              %%Keyword: Thresholds <advanced>
              <HELP>
              Specifies convergence thresholds. Four individual thresholds are specified
              as arguments, which have to be fulfilled simultaneously to reach convergence:
              EThr, DThr and FThr
              specify the maximum permissible difference in energy, density matrix elements
              and Fock matrix elements, respectively, in the last two iterations. The
              DltNTh finally specifies the norm of the orbital displacement vector used
              for the orbital rotations in the second-order/C2-DIIS procedure.
              The corresponding values are read in the order given above.
              The default values are 1.0d-9, 1.0d-4, 1.5d-4, and 0.2d-4,
              respectively.
              Note that these thresholds automatically define the threshold
              used in the direct Fock matrix construction to estimate individual
              contributions to the Fock matrix such that
              the computed energy will have an accuracy that is better than the
              convergence threshold.
              </HELP>
              </KEYWORD>

:kword:`NODIis`
  Disable the DIIS convergence acceleration procedure.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="NODIIS" APPEAR="NoDIIS" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: NoDIIS <advanced>
              <HELP>
              Disable the DIIS convergence acceleration procedure.
              </HELP>
              </KEYWORD>

:kword:`DIISthr`
  Set the threshold on the change in density, at which the DIIS procedure
  is turned on.
  The keyword takes as argument a (double precision) floating point number.
  The default value is 0.15.

  .. xmldoc:: %%Keyword: DIISthr <advanced>
              Set the threshold on the change in density, at which the DIIS procedure
              is turned on.
              The keyword takes as argument a (double precision) floating point number.
              The default value is 0.15.

:kword:`QNRThr`
  Set the threshold on the change in density, at which the
  second-order/\ :math:`C^2`\-DIIS
  procedure kicks in.
  The keyword takes as argument a (double precision) floating point number.
  The default value is 0.15.

  **Note:** the change in density has to drop under both the
  :kword:`DIISthr` and the :kword:`QNRThr` threshold, for the
  second-order/\ :math:`C^2`\-DIIS to be activated. If the latter is set to zero
  the older first order :math:`C^2`\-DIIS procedure will be used instead.

  .. xmldoc:: %%Keyword: QNRTHR <advanced>
              Set the threshold on the change in density, at which the
              second-order/C2-DIIS
              procedure kicks in.
              The keyword takes as argument a (double precision) floating point number.
              The default value is 0.15.
              Note: the change in density has to drop under both the
              DIISthr and the QNRThr threshold, for the
              second-order/C2-DIIS to be activated. If the latter is set to zero
              the older first order C2-DIIS procedure will be used instead.

:kword:`C1DIis`
  Use :math:`C^1`\-DIIS for convergence acceleration rather than :math:`C^2`\-DIIS
  which is the default (not recommended).

  .. xmldoc:: %%Keyword: C1DIIS <advanced>
              Use C1-DIIS for convergence acceleration rather than C2-DIIS
              which is the default (not recommended).

:kword:`NODAmp`
  Disable the Damping convergence acceleration procedure.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="NODAMP" APPEAR="NoDAMP" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: NODAMP <advanced>
              <HELP>
              Disable the Damping convergence acceleration procedure.
              </HELP>
              </KEYWORD>

:kword:`OCCNumbers`
  Gives the option to specify occupation numbers other than 0 and 2.
  This can be useful for generating starting orbitals for open shell
  cases. It should be noted however, that it is still the closed shell
  SCF energy functional that is optimized, thus yielding unphysical
  energies. Occupation numbers have to be provided for all occupied
  orbitals.
  In the case of UHF calculation occupation numbers should be specified
  on two different entries: for alpha and beta spin.

  .. xmldoc:: %%Keyword: Occnumbers <advanced>
              Gives the option to specify occupation numbers other than 0 and 2.
              This can be useful for generating starting orbitals for open shell
              cases. It should be noted however, that it is still the closed shell
              SCF energy functional that is optimized, thus yielding unphysical
              energies. Occupation numbers have to be provided for all occupied
              orbitals.
              In the case of UHF calculation occupation numbers should be specified
              on two different entries: for alpha and beta spin

:kword:`IVO`
  Specifies that the virtual orbitals are to be improved for
  subsequent MCSCF calculations. The core Hamiltonian is diagonalized
  within the virtual orbital subspace, thus yielding as compact orbitals
  as possible with the constraint that they have to be orthogonal to the
  occupied orbitals.
  **Note** that this option must not be used whenever the Hartree--Fock
  wavefunction itself is used as a reference in a subsequent calculation.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="IVO" APPEAR="IVO" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: IVO <advanced> GUI:keyword
              <HELP>
              Specifies that the virtual orbitals are to be improved for
              subsequent MCSCF calculations. The core Hamiltonian is diagonalized
              within the virtual orbital subspace, thus yielding as compact orbitals
              as possible with the constraint that they have to be orthogonal to the
              occupied orbitals.
              Note that this option must not be used whenever the Hartree-Fock
              wavefunction itself is used as a reference in a subsequent calculation.
              </HELP>
              </KEYWORD>

:kword:`NOMInimization`
  Program will use density differences
  :math:`D^{(k)}-D^{(k-1)}`
  rather than minimized differences.

  .. xmldoc:: %%Keyword: nominimization <advanced>
              Program will use density differences
              D(k)-D(k-1)
              rather than minimized differences.

:kword:`ONEGrid`
  Disable use of a smaller intermediate grid in the integration of the
  DFT functional during the first SCF iterations.

  .. xmldoc:: <KEYWORD MODULE="SCF" NAME="ONEGRID" APPEAR="OneGrid" KIND="SINGLE" LEVEL="ADVANCED" REQUIRE="KSDFT">
              %%Keyword: ONEGrid <advanced>
              <HELP>
              Disable use of a smaller intermediate grid in the integration of the
              DFT functional during the first SCF iterations.
              </HELP>
              </KEYWORD>

:kword:`RFPErt`
  This keyword will add a constant reaction field perturbation to the
  bare nuclei hamiltonian.
  The perturbation is read from :file:`RUNOLD` (if not present defaults to :file:`RUNFILE`) and
  is the latest self consistent perturbation generated
  by one of the programs :program:`SCF` or :program:`RASSCF`.

  .. xmldoc:: %%Keyword: Rfpert <advanced>
              This keyword will add a constant reaction field perturbation to the
              bare nuclei hamiltonian.
              The perturbation is read from RUNOLD (if not present defaults to RUNFILE) and
              is the latest selfconsistent perturbation generated
              by one of the programs SCF or RASSCF.

:kword:`STAT`
  This keyword will add an addition print outs with statistic information.

  .. xmldoc:: %%Keyword: STAT <advanced>
              This keyword will add an addition print outs with statistic information

For calculations of a molecule in a reaction field see section :ref:`UG:sec:rfield`
of the present manual and section :ref:`TUT:sec:cavity` of the examples manual.

.. include:: ../dft_functionals.inc

Keywords for direct calculations
................................

*Note* again that the threshold for contributions to
the Fock matrix depends on the convergence thresholds
mentioned above. The choice between the conventional and direct SCF
methods is based on the presence of a two-electron integral file
(file :file:`ORDINT`). The keyword
:kword:`Direct` in the :program:`SEWARD` input controls that no
two-electron integral file is to be generated and that integral direct
algorithms can be used in subsequent modules. Thus, *the
choice between conventional and direct SCF is done already in the input
for the integral program* :program:`SEWARD`. The direct (or semi-direct)
path will be taken whenever there are no two-electron integrals available.

.. class:: keywordlist

:kword:`CONVentional`
  This option will override the automatic choice between the conventional
  and the direct SCF algorithm such that the conventional method will
  be executed regardless of the status of the :file:`ORDINT` file.

  .. xmldoc:: %%Keyword: conventional <advanced>
              This option will override the automatic choice between the conventional
              and the direct SCF algorithm such that the conventional method will
              be executed regardless of the status of the ORDINT file.

:kword:`DISK`
  This option enables/disables the semi-direct algorithm. It requires
  two arguments which specifies the max Mbyte of integrals that are written
  on disk during the first iteration (and retrieved later in subsequent
  iterations)
  and the size of the corresponding I/O buffer in kbyte.
  The default values are
  2000 MByte and 512 kByte. In case the specified disk space is zero and the I/O buffer
  is different from zero it will default to a semi-direct SCF with in-core storage
  of the integrals. The size of the memory for integrals storage is the size of the
  I/O buffer. If the size of the disk is non-zero and the I/O buffer size is zero the
  latter will be reset to the default value.

  .. xmldoc:: %%Keyword: Disk <advanced>
              This option enables/disables the semi-direct algorithm. It requires
              two arguments which specifies the max Mbyte of integrals that are written
              on disk during the first iteration (and retrieved later in subsequent
              iterations) and the size of the I/O buffer in kbyte.
              The default values are
              2000 MByte and 512 kByte. In case the specified disk space is zero and the I/O buffer
              is different from zero it will default to a semi-direct SCF with in-core storage
              of the integrals. The size of the memory for integrals storage is the size of the
              corresponding I/O buffer. If the size of the disk is non-zero and the I/O buffer size is zero the
              latter will be reset to the default value.

:kword:`THIZe`
  This option specifies a threshold for two-electron integrals.
  Only integrals above this threshold (but not necessarily all of those) are kept
  on disk for the semi-direct algorithm.
  The keyword takes as argument a (double precision) floating point number.

  .. xmldoc:: %%Keyword: Thize <advanced>
              This option specifies a threshold for two-electron integrals.
              Only integrals above this threshold (but not necessarily all of those) are kept
              on disk for the semi-direct algorithm.
              The keyword takes as argument a (double precision) floating point number.

:kword:`SIMPle`
  If this option is specified, only a simple prescreening scheme,
  based solely on the
  estimated two-electron integral value will be employed (no density involved).

  .. xmldoc:: %%Keyword: Simple <advanced>
              If this option is specified, only a simple prescreening scheme,
              based solely on the
              estimated two-electron integral value will be employed (no density involved).

Limitations
...........

The limitations/MODULE on the number of basis functions are the same as specified
for :program:`SEWARD`.

Input examples
..............

First we have the bare minimum of input. This will work well for almost
all systems containing an even number of electrons. ::

  &SCF

The next example is almost as simple. Here we have an open shell case,
i.e. you have an odd number of electrons in the neutral system and you
need to generate starting orbitals for :program:`RASSCF`.
In this case we recommend that you perform a calculation on the
cation with the input below. ::

  &SCF; Charge= 1

The next example explains how to run UHF code for a nitrogen atom: ::

  &SCF; UHF; ZSPIN=3

The next example is a bit more elaborate and show how to use
a few of the keywords. The system is water that have the
electron configuration :math:`\text{1a}_1^2 \text{2a}_1^2 \text{3a}_1^2 \text{1b}_1^2 \text{1b}_2^2`. ::

  &SCF; Title= Water molecule. Experimental equilibrium geometry. The symmetries are a1, b2, b1 and a2.
  Occupied= 3 1 1 0
  Threshold= 0.5D-9 0.5D-6 0.5D-6 0.5D-5
  * semi-direct algorithm writing max 128k words (1MByte) to disk
  * the size of the I/O buffer by default (512 kByte)
  Disk= 1 0
  Ivo

.. xmldoc:: </MODULE>
