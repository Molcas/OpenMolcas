.. index::
   single: Program; CASVB
   single: CASVB

.. _UG\:sec\:casvb:

:program:`casvb`
================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. The input format of CASVB is not completely compatible with the XML specification
   format used by MolGUI.

.. xmldoc:: <MODULE NAME="CASVB">
            %%Description:
            <HELP>
            CASVB is a program for performing general valence bond calculations.
            It can be used in two basic modes:
            a) variational optimization of quite general types of
               nonorthogonal MCSCF or modern valence bond wavefunctions, or
            b) representation of CASSCF wavefunctions in modern valence form,
               using overlap- (relatively inexpensive) or energy-based criteria.
            </HELP>

This program can be used in two basic modes:

#. variational optimization of quite general types of
   nonorthogonal MCSCF or modern valence bond wavefunctions
#. representation of CASSCF wavefunctions in modern valence form,
   using overlap- (*relatively inexpensive*) or energy-based criteria.

For generating representations of CASSCF wavefunctions, the program
is invoked by the command :kword:`CASVB`.
For variational optimization of wavefunctions it is normally invoked
inside :program:`RASSCF` by the sub-command :kword:`VB` (see :ref:`the RASSCF documentation <vbinrasscf>`).

Bibliography: see :cite:`casvb1,casvb2,casvb3,casvb4`.

.. index::
   pair: Dependencies; CASVB

.. _UG\:sec\:casvb_dependencies:

Dependencies
------------

The :program:`CASVB` program needs the :file:`JOBIPH` file from a :program:`RASSCF` calculation,
and in addition also the :file:`ONEINT` and :file:`ORDINT` files from :program:`SEWARD`.

.. index::
   pair: Files; CASVB

.. _UG\:sec\:casvb_files:

Files
-----

Input files
...........

:program:`CASVB` will use the following input
files: :file:`ONEINT`, :file:`ORDINT`, :file:`RUNFILE`, :file:`JOBIPH`,
(for more information see :numref:`UG:sec:files_list`), and
:file:`VBWFN` with
valence bond wavefunction information (orbital and structure coefficients).

Output files
............

.. class:: filelist

:file:`JOBIPH`
  On exit, the :program:`RASSCF` interface file is overwritten with the
  CASVB wavefunction.

:file:`VBWFN`
  Valence bond wavefunction information (orbital and structure coefficients).

.. _UG\:sec\:casvb_input:

.. index::
   pair: Input; CASVB

Input
-----

This section describes the input to the :program:`CASVB` program.
The input for each module is preceded by its name like: ::

 &CASVB

.. index::
   pair: Keywords; CASVB

Keywords
........

Optional keywords

.. class:: keywordlist

:kword:`END of Input`
  This marks the end of the input to the program.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="END" APPEAR="End of input" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: END of Input <basic>
              <HELP>
              This marks the end of the input to the program and is the only compulsory
              keyword.
              </HELP>
              </KEYWORD>

Optional keywords to define the CASSCF wavefunction. Not generally required
because values stored in the job interface
file or used by the :program:`RASSCF` program will normally be appropriate.

.. class:: keywordlist

:kword:`FROZen`
  Specifies frozen orbitals, as in the :program:`RASSCF` program.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="FROZEN" APPEAR="Frozen orbitals" KIND="INTS_LOOKUP" SIZE="NSYM" LEVEL="ADVANCED" MIN_VALUE="0">
              %%Keyword: FROZen <advanced>
              <HELP>
              Specifies frozen orbitals, as in the RASSCF program. This
              keyword is generally not required because the value stored in the job
              interface file or used by the RASSCF program will normally be appropriate.
              </HELP>
              </KEYWORD>

:kword:`INACtive`
  Specifies inactive orbitals, as in the :program:`RASSCF` program.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="INACTIVE" APPEAR="Inactive orbitals" KIND="INTS_LOOKUP" SIZE="NSYM" LEVEL="ADVANCED" MIN_VALUE="0">
              %%Keyword: INACtive <advanced>
              <HELP>
              Specifies inactive orbitals, as in the RASSCF program. This
              keyword is generally not required because the value stored in the job
              interface file or used by the RASSCF program will normally be appropriate.
              </HELP>
              </KEYWORD>

:kword:`NACTel`
  Specifies the number of active electrons, as in the :program:`RASSCF` program.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="NACTEL" APPEAR="Active electrons" KIND="INTS" SIZE="3" LEVEL="ADVANCED" MIN_VALUE="0">
              %%Keyword: NACTel <advanced>
              <HELP>
              Specifies number of active electrons, as in the RASSCF program. This
              keyword is generally not required because the value stored in the job
              interface file or used by the RASSCF program will normally be appropriate.
              </HELP>
              </KEYWORD>

:kword:`RAS2`
  Specifies RAS2 orbitals, as in the :program:`RASSCF` program.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="RAS2" APPEAR="RAS2" LEVEL="ADVANCED" KIND="INTS_LOOKUP" SIZE="NSYM" DEFAULT_VALUE="0" MIN_VALUE="0">
              %%Keyword: RAS2 <advanced>
              <HELP>
              Specifies RAS2 orbitals, as in the RASSCF program. This
              keyword is generally not required because the value stored in the job
              interface file or used by the RASSCF program will normally be appropriate.
              </HELP>
              </KEYWORD>

:kword:`SPIN`
  Specifies the total spin, as in the :program:`RASSCF` program.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="SPIN" APPEAR="Spin" LEVEL="ADVANCED" KIND="INT" DEFAULT_VALUE="1" MIN_VALUE="1">
              %%Keyword: SPIN <advanced>
              <HELP>
              Specifies the total spin, as in the RASSCF program. This
              keyword is generally not required because the value stored in the job
              interface file or used by the RASSCF program will normally be appropriate.
              </HELP>
              </KEYWORD>

:kword:`SYMMetry`
  Specifies the CASSCF wavefunction symmetry, as in the :program:`RASSCF` program.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="SYMMETRY" APPEAR="Symmetry" LEVEL="ADVANCED" KIND="INT" DEFAULT_VALUE="1" MIN_VALUE="1" MAX_VALUE="8">
              %%Keyword: SYMMetry <advanced>
              <HELP>
              Specifies the CASSCF wavefunction symmetry, as in the RASSCF program. This
              keyword is generally not required because the value stored in the job
              interface file or used by the RASSCF program will normally be appropriate.
              </HELP>
              </KEYWORD>

Optional keywords to define the VB wavefunction

.. class:: keywordlist

:kword:`CON`
  .. index::
     single: CON

  The spatial
  VB configurations are defined in terms of the active orbitals, and may be
  specified using one or more :kword:`CON` keywords: ::

    CON
    n1 n2 n3 n4 ...

  The configurations can be specified by occupation numbers, so that
  :math:`n_i` is the occupation of the :math:`i`\th valence bond orbital. Alternatively a list of
  :math:`N_{\text{act}}` orbital numbers (in any order) may be provided --- the
  program determines which definition applies. The two specifications ``1 0 1 2``
  and ``1 3 4 4`` are thus equivalent.

  Input configurations are reordered by :program:`CASVB`, so that configurations have
  non-decreasing double occupancies. Configurations that are inconsistent with the
  value for the total spin are ignored.

  If no configurations are specified the single "covalent" configuration
  :math:`\phi_1\phi_2\cdots\phi_{N_{\text{act}}}` is assumed.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="CON" APPEAR="Configurations" LEVEL="BASIC" KIND="STRING">
              %%Keyword: CON <basic>
              <HELP>
              Specifies spatial VB configurations in terms of the active orbitals.
              The default is a single configuration of singly-occupied orbitals.
              </HELP>
              </KEYWORD>

:kword:`COUPle`
  .. index::
     single: COUPLE

  ::

    COUPLE
    key

  ``key`` may be chosen from ``KOTANI`` (default), ``SERBER``, ``RUMER``, ``PROJECT`` or ``LTRUMER``,
  specifying the scheme for constructing the
  spin eigenfunctions used in the definition of valence bond structures. ``PROJECT``
  refers to spin functions generated using a spin projection operator, ``LTRUMER`` to
  Rumer functions with the so-called "leading term" phase convention.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="COUPLE" APPEAR="Couple scheme" LEVEL="ADVANCED" KIND="CHOICE" LIST="KOTANI,SERBER,RUMER,PROJECT,LTRUMER" DEFAULT_VALUE="KOTANI">
              %%Keyword: COUPle <advanced>
              <HELP>
              Specifies the scheme for constructing the spin eigenfunctions to be used.
              Possible values: KOTANI (default), SERBER, RUMER, PROJECT, LTRUMER
              </HELP>
              </KEYWORD>

:kword:`WAVE`
  .. index::
     single: WAVE

  ::

    WAVE
    N S1 S2 ...

  This keyword can be used to specify explicitly the number of electrons and spin(s) to
  be used with a configuration list. If :math:`N` is less than the present number of active electrons,
  the input wavefunction fragment is assumed to form part of a direct product. Otherwise, the spins
  specified may be greater than or equal to the :kword:`SPIN` value specified as input to the :program:`RASSCF`
  program. Defaults, for both :math:`N` and :math:`S`, are the values used by :program:`RASSCF`.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="WAVE" APPEAR="Wavefunction" LEVEL="ADVANCED" KIND="STRING">
              %%Keyword: WAVE <advanced>
              <HELP>
              Specifies number of electrons and spins to be used with a configuration list.
              Defaults are the values used by RASSCF.
              </HELP>
              </KEYWORD>

Optional keywords for the recovery and/or storage of orbitals and vectors

.. class:: keywordlist

:kword:`STARt`
  .. index::
     single: START

  ::

    START
    key-1=filename-1
    key-2=filename-2
    ...

  Specifies input files for VB wavefunction (``key-i``\=VB),
  CASSCF CI vector (``key-i``\=CI) and/or CASSCF molecular orbitals
  (``key-i``\=MO).
  By default, the required information is taken from the file :file:`JOBOLD`.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="START" APPEAR="Input files" LEVEL="ADVANCED" KIND="UNKNOWN">
              %%Keyword: STARt <advanced>
              <HELP>
              Specifies various input files. Default is to take the required information
              from JOBOLD.
              </HELP>
              </KEYWORD>

:kword:`SAVE`
  .. index::
     single: SAVE

  ::

    SAVE
    key-1=filename-1
    key-2=filename-2
    ...

  Specifies output files for VB wavefunction (``key-i``\=VB)
  and/or the VB CI vector (``key-i``\=VBCI). By default, the VB CI
  vector is written to the file JOBIPH.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="SAVE" APPEAR="Output files" LEVEL="ADVANCED" KIND="UNKNOWN">
              %%Keyword: SAVE <advanced>
              <HELP>
              Specifies output files. By default, the VB CI vector is written to the
              file JOBIPH.
              </HELP>
              </KEYWORD>

Optional keywords to override the starting guess

.. class:: keywordlist

:kword:`GUESs`
  .. index::
     single: GUESS

  ::

    GUESS
    key-1 ...
    key-2 ...
    ENDGUESs

  The :kword:`GUESS` keyword initiates the input of a guess for the valence bond orbitals and/or
  structure coefficients. ``key-i`` can be either :kword:`ORB` or :kword:`STRUC`.
  These keywords
  modify the guess provided by the program. It is
  thus possible to modify individual orbitals in a previous solution
  so as to construct the starting
  guess. The :kword:`ENDGUESs` keyword terminates the guess input. ::

    ORB
    i c1 c2 ... cmact

  Specifies a starting guess for valence bond orbital number :math:`i`. The guess is specified
  in terms of the :math:`m_{\text{act}}` active MOs defining the CASSCF wavefunction. ::

    STRUC
    c1 c2 ... cNVB

  Specifies a starting guess for the :math:`N_{\text{VB}}` structure coefficients. If this keyword
  is not provided, the perfect-pairing mode of
  spin coupling is assumed for the spatial configuration having the least
  number of doubly occupied orbitals.
  Note that the definition of structures depends on the value of :kword:`COUPLE`. Doubly occupied
  orbitals occur first in all configurations, and the spin eigenfunctions are based on the singly
  occupied orbitals being in ascending order.

  .. xmldoc:: <GROUP MODULE="CASVB" NAME="GUESS" APPEAR="Guess" KIND="BLOCK" LEVEL="BASIC">

  .. xmldoc:: %%Keyword: GUESs <basic>
              <HELP>
              Initiates guess input. Sub-keywords are ORB and STRUC, as described in the
              manual. The ENDGUESs keyword terminates the guess input.
              </HELP>

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="ORB" APPEAR="Orbital" LEVEL="BASIC" KIND="CUSTOM" />

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="STRUC" APPEAR="Structure" LEVEL="BASIC" KIND="UNKNOWN" />

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="AOBASIS" LEVEL="UNDOCUMENTED" KIND="SINGLE" />

  .. xmldoc:: </GROUP>

:kword:`ORBPerm`
  .. index::
     single: ORBPERM

  ::

    ORBPERM
    i1 ... imact

  Permutes the orbitals in the valence bond wavefunction and changes their phases according to
  :math:`\phi_j'=\sign(i_j)\phi_{\abs(i_j)}`. The guess may be further modified using the
  :kword:`GUESS` keyword. Additionally, the structure coefficients will be transformed
  according to the given permutation (note that the configuration list must be closed under
  the orbital permutation for this to be possible).

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="ORBPERM" APPEAR="Orbital permutation" LEVEL="ADVANCED" KIND="STRING">
              %%Keyword: ORBPerm <advanced>
              <HELP>
              Permutes the VB orbitals, and modifies phases, as described in the manual.
              </HELP>
              </KEYWORD>

Optional keywords for optimization control

.. class:: keywordlist

:kword:`CRIT`
  .. index::
     single: CRIT

  ::

    CRIT
    method

  .. compound::

    Specifies the criterion for the optimization. ``method`` can be :kword:`OVERLAP` or :kword:`ENERGY`
    (:kword:`OVERLAP` is default).
    The former maximizes the normalized overlap with the CASSCF wavefunction:

    .. math:: \max\left(\frac{\braket{\Psi_{\text{CAS}}}{\Psi_{\text{VB}}}} {\left(\braket{\Psi_{\text{VB}}}{\Psi_{\text{VB}}}\right)^{1/2}}\right)

    and the latter simply minimizes the energy:

    .. math:: \min\left(\frac{\braopket{\Psi_{\text{VB}}}{\hat{H}}{\Psi_{\text{VB}}}}{\braket{\Psi_{\text{VB}}}{\Psi_{\text{VB}}}}\right).

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="CRIT" APPEAR="Optimization criterion" LEVEL="BASIC" KIND="CHOICE" LIST="OVERLAP,ENERGY" DEFAULT_VALUE="OVERLAP">
              %%Keyword: CRIT <basic>
              <HELP>
              Defines the optimization criterion.
              Possible values: OVERLAP (default) or ENERGY.
              </HELP>
              </KEYWORD>

:kword:`MAXIter`
  .. index::
     single: MAXITER

  ::

    MAXITER
    Niter

  Specifies the maximum number of iterations in the second-order optimizations. Default is :math:`N_{\text{iter}}`\=50.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="MAXITER" APPEAR="Maximum iterations" LEVEL="ADVANCED" KIND="INT" DEFAULT_VALUE="50" MIN_VALUE="1">
              %%Keyword: MAXIter <advanced>
              <HELP>
              Specifies the maximum number of iterations to be used. Default value is 50.
              </HELP>
              </KEYWORD>

:kword:`(NO)CASProj`
  .. index::
     single: CASPROJ
     single: NOCASPROJ

  ::

    (NO)CASPROJ

  With this keyword the structure coefficients are picked from the transformed CASSCF CI vector, leaving
  only the orbital variational parameters. For further details see the bibliography.
  This option may be useful to aid convergence.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="CASPROJ" APPEAR="CAS proj" LEVEL="ADVANCED" KIND="SINGLE" EXCLUSIVE="NOCASPROJ,NOPROJCAS">
              %%Keyword: CASProj <advanced>
              <HELP>
              Defines structure coefficients from transformed CASSCF wavefunction.
              </HELP>
              </KEYWORD>

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="NOCASPROJ" APPEAR="No CAS proj" LEVEL="ADVANCED" KIND="SINGLE" EXCLUSIVE="CASPROJ,PROJCAS">
              %%Keyword: NOCASProj <advanced>
              <HELP>
              Disables CASProj
              </HELP>
              </KEYWORD>

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="PROJCAS" LEVEL="UNDOCUMENTED" KIND="SINGLE" EXCLUSIVE="NOCASPROJ,NOPROJCAS" />

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="NOPROJCAS" LEVEL="UNDOCUMENTED" KIND="SINGLE" EXCLUSIVE="CASPROJ,PROJCAS" />

:kword:`SADDle`
  .. index::
     single: SADDLE

  ::

    SADDLE
    n

  Defines optimization onto an :math:`n`\th-order saddle point.
  See also :cite:`casvb5`.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="SADDLE" APPEAR="Saddle point" LEVEL="ADVANCED" KIND="INT" MIN_VALUE="1">
              %%Keyword: SADDLe <advanced>
              <HELP>
              Defines optimization onto an n-th order saddle point.
              </HELP>
              </KEYWORD>

:kword:`(NO)INIT`
  .. index::
     single: INIT
     single: NOINIT

  ::

    (NO)INIT`

  Requests a sequence of preliminary optimizations which aim to minimize the
  computational cost while maximizing the likelihood of stable
  convergence. This feature is the default if no wavefunction guess is available
  and no :kword:`OPTIM` keyword specified in the input.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="INIT" APPEAR="Initial optimizations" LEVEL="ADVANCED" KIND="SINGLE" EXCLUSIVE="NOINIT">
              %%Keyword: INIT <advanced>
              <HELP>
              Requests a sequence of preliminary optimizations which aim to minimize
              the computational cost while maximizing the likelihood of stable
              convergence. This is the default behaviour when no wavefunction guess is
              available and no OPTIM keyword has been specified.
              </HELP>
              </KEYWORD>

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="NOINIT" APPEAR="No initial optimizations" LEVEL="ADVANCED" KIND="SINGLE" EXCLUSIVE="INIT">
              %%Keyword: NOINIT <advanced>
              <HELP>
              Disables INIT
              </HELP>
              </KEYWORD>

:kword:`METHod`
  .. index::
     single: METHOD

  ::

    METHOD
    key

  Selects the optimization algorithm to be used. ``key`` can be one
  of: :kword:`FLETCHER`, :kword:`TRIM`, :kword:`TRUSTOPT`, :kword:`DAVIDSON`,
  :kword:`STEEP`, :kword:`VB2CAS`, :kword:`AUGHESS`, :kword:`AUG2`,
  :kword:`CHECK`, :kword:`DFLETCH`, :kword:`NONE`, or :kword:`SUPER`. Recommended are
  the direct procedures :kword:`DFLETCH` or :kword:`AUGHESS`. For general
  saddle-point optimization :kword:`TRIM` is used. Linear (CI only) optimization
  problems use :kword:`DAVIDSON`. :kword:`NONE` suspends optimization, while
  :kword:`CHECK` carries out a finite-difference check of the gradient and Hessian.

  The default algorithm chosen by :program:`CASVB` will be usually be adequate.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="METHOD" APPEAR="Method" LEVEL="ADVANCED" KIND="CHOICE" LIST="FLETCHER,TRIM,TRUSTOPT,DAVIDSON,STEEP,VB2CAS,AUGHESS,AUG2,CHECK,DFLETCH,NONE,SUPER">
              %%Keyword: METHod <advanced>
              <HELP>
              Selects optimization algorithm.
              Possible values: FLETCHER, TRIM, TRUSTOPT, DAVIDSON, STEEP, VB2CAS,
              AUGHESS, AUG2, CHECK, DFLETCH, NONE or SUPER.
              The default algorithm chosen by CASVB will usually be adequate.
              </HELP>
              </KEYWORD>

:kword:`TUNE`
  .. index::
     single: TUNE

  ::

    TUNE
    ...

  Enables the input of individual parameters to be used in the optimization procedure
  (*e.g.* for controlling step-size selection and convergence testing).
  Details of the values used are output if ``print(3)`` :math:`\geq` 3 is specified.
  For expert use only.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="TUNE" APPEAR="Tune" LEVEL="ADVANCED" KIND="UNKNOWN">
              %%Keyword: TUNE <advanced>
              <HELP>
              Enables the input of individual parameters to be used in the optimization procedure. Expert use only. See manual.
              </HELP>
              </KEYWORD>

:kword:`OPTIm`
  .. index::
     single: OPTIM

  More than one optimization may be performed in the same :program:`CASVB` run,
  by the use of :kword:`OPTIM` keywords: ::

    OPTIM
    [...
    ENDOPTIM]

  The subcommands may be any optimization declarations defined in this
  section, as well as any symmetry or constraints specifications.
  Commands given as arguments to :kword:`OPTIM`
  will apply only to this optimization step, whereas commands specified
  outside will act as default definitions for all subsequent :kword:`OPTIM`
  specifications.

  The :kword:`OPTIM` keyword
  need not be specified if only one optimization step is required,

  When only a machine-generated guess is available, :program:`CASVB` will
  attempt to
  define a sequence of optimization steps that aims to maximize the
  likelihood of successful convergence (while minimizing
  CPU usage). To override this behaviour, simply specify one or more
  :kword:`OPTIM` keywords. The :kword:`ENDOPTIm` keyword marks the end of the
  specifications of an optimization step.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="OPTIM" APPEAR="Optimizations" LEVEL="ADVANCED" KIND="SINGLE">
              %%Keyword: OPTIm <advanced>
              <HELP>
              Defines one or more optimization steps. Subcommands can be any
              optimization declarations, as well as any symmetry or constraints
              specifications. Usually omitted if only one optimization step is required.
              Terminated by the keyword ENDOPTIm.
              </HELP>
              </KEYWORD>

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="ENDOPTIM" LEVEL="UNDOCUMENTED" KIND="SINGLE" />

:kword:`ALTErn`
  .. index::
     single: ALTERN

  A loop over two or more optimization steps may be specified using: ::

    ALTERN
    Niter
    ...
    ENDALTERN

  The program will repeat the specified optimization steps
  until either all optimizations have converged, or the maximum iteration count,
  :math:`N_{\text{iter}}`, has been reached.
  The :kword:`ENDALTErn` keyword marks the end of the specification of an
  ALTERN loop.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="ALTERN" APPEAR="Alternate" LEVEL="ADVANCED" KIND="UNKNOWN">
              %%Keyword: ALTErn <advanced>
              <HELP>
              Defines alternating optimizations over two or more optimization steps (see
              manual). Terminated by the ENDALTErn keyword.
              </HELP>
              </KEYWORD>

Optional keywords for definitions of molecular symmetry and any
constraints on the VB wavefunction

.. class:: keywordlist

:kword:`SYMElm`
  .. index::
     single: SIMELM

  Various issues associated with symmetry-adapting valence bond wavefunctions
  are considered, for example, in :cite:`casvb6`. ::

    SYMELM
    label sign

  Initiates the definition of a symmetry operation referred to by ``label`` (any three characters).
  ``sign`` can be :math:`+` or :math:`-`; it specifies whether the total wavefunction is symmetric or
  antisymmetric under this operation, respectively. A value for ``sign`` is not always necessary
  but, if provided, constraints will be put on the structure coefficients to ensure that the
  wavefunction has the correct overall symmetry (note that the configuration list must be closed
  under the orbital permutation induced by ``label`` for this to be possible).
  The default for ``label`` is the identity.

  The operator is defined in terms of its action on the active MOs as specified by
  one or more of the keywords :kword:`IRREPS`, :kword:`COEFFS`, or :kword:`TRANS`. Any
  other keyword, including optional use of the :kword:`ENDSYMElm` keyword, will
  terminate the definition of this symmetry operator. ::

    IRREPS
    i1 i2 ...

  The list :math:`i_1, i_2 \ldots` specifies which irreducible representations (as defined in
  the CASSCF wavefunction) are antisymmetric with respect to the ``label`` operation.
  If an irreducible representation is not otherwise specified it is assumed to be symmetric
  under the symmetry operation. ::

    COEFFS
    i1 i2 ...

  The list :math:`i_1, i_2 \ldots` specifies which individual CASSCF MOs are antisymmetric with
  respect to the ``label`` operation. If an MO is not otherwise specified, it is assumed to be
  symmetric under the symmetry operation. This specification may be useful if, for example, the
  molecule possesses symmetry higher than that exploited in the CASSCF calculation. ::

    TRANS
    ndim i1 ... indim c1,1 c1,2 ... cndim,ndim

  Specifies a general :math:`n_{\text{dim}}\times n_{\text{dim}}` transformation involving the MOs :math:`i_1,
  \ldots i_{n_{\text{dim}}}`,
  specified by the :math:`c` coefficients. This may be useful for systems with a two- or
  three-dimensional irreducible representation, or if localized orbitals define the CASSCF
  wavefunction. Note that the specified transformation must always be orthogonal.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="SYMELM" APPEAR="Symmetry elements" LEVEL="ADVANCED" KIND="STRING">
              %%Keyword: SYMElm <advanced>
              <HELP>
              Initiates the definition of a symmetry operation (see manual).
              Sub-keywords are IRREPS, COEFFS, or TRANS. Terminated with ENDSYMElm.
              </HELP>
              </KEYWORD>

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="IRREPS" LEVEL="UNDOCUMENTED" KIND="INTS_LOOKUP" SIZE="ANY" />

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="COEFFS" LEVEL="UNDOCUMENTED" KIND="INTS_LOOKUP" SIZE="ANY" />

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="ENDSYMELM" LEVEL="UNDOCUMENTED" KIND="SINGLE" />

:kword:`ORBRel`
  .. index::
     single: ORBREL

  In general, for a VB wavefunction to be symmetry-pure, the orbitals must form a representation
  (not necessarily irreducible) of the symmetry group. Relations between orbitals under
  the symmetry operations defined by :kword:`SYMELM` may be specified according to: ::

    ORBREL
    i1 i2 label-1 label-2 ...

  Orbital :math:`i_1` is related to orbital :math:`i_2` by the sequence of operations defined by the ``label``
  specifications (defined previously using :kword:`SYMELM`). The operators operate right to left. Note
  that :math:`i_1` and :math:`i_2` may coincide. Only the minimum number of
  relations required to define all the orbitals should be provided; an error exit
  will occur if redundant :kword:`ORBREL` specifications are found.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="ORBREL" APPEAR="Orbital relations" LEVEL="ADVANCED" KIND="CUSTOM">
              %%Keyword: ORBRel <advanced>
              <HELP>
              Specifies the relationship between two VB orbitals under symmetry
              operation(s) defined by SYMElm. See manual.
              </HELP>
              </KEYWORD>

:kword:`(NO)SYMProj`
  .. index::
     single: SYMPROJ
     single: NOSYMPROJ

  As an alternative to incorporating constraints, one may also ensure correct
  symmetry of the wavefunction by use of a projection operator: ::

    (NO)SYMPROJ
    [irrep-1 irrep-2 ...]

  The effect of this keyword is to set to zero the coefficients in unwanted
  irreducible representations.
  For this purpose, the symmetry group defined for the CASSCF wavefunction
  is used (always a subgroup of :math:`D_{2h}`).
  The list of irreps in the command specifies which components
  of the wavefunction should be kept.
  If no irreducible representations are given, the current
  wavefunction symmetry is assumed. In a state-averaged calculation,
  all irreps are retained for which a non-zero weight has been specified in the
  wavefunction definition.
  The :kword:`SYMPROJ` keyword may also be used in combination with constraints.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="SYMPROJ" APPEAR="Symmetry projection" LEVEL="ADVANCED" KIND="INTS_LOOKUP" SIZE="NSYM" EXCLUSIVE="NOSYMPROJ">
              %%Keyword: SYMProj <advanced>
              <HELP>
              Projects the VB wavefunction onto given irrep(s). See manual.
              </HELP>
              </KEYWORD>

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="NOSYMPROJ" APPEAR="No symmetry projection" LEVEL="ADVANCED" KIND="INTS_LOOKUP" SIZE="NSYM" EXCLUSIVE="SYMPROJ">
              %%Keyword: NoSYMProj <advanced>
              <HELP>
              Disables SYMProj
              </HELP>
              </KEYWORD>

:kword:`FIXOrb`
  .. index::
     single: FIXORB

  ::

    FIXORB
    i1 i2 ...

  This command freezes the orbitals specified in the list
  :math:`i_1, i_2 \ldots` to that of the starting guess. Alternatively the
  special keywords :kword:`ALL` or :kword:`NONE` may be used. These orbitals
  are eliminated from the optimization procedure, but will still be
  normalized and symmetry-adapted according to any :kword:`ORBREL`
  keywords given.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="FIXORB" APPEAR="Freeze orbitals" LEVEL="ADVANCED" KIND="INTS_LOOKUP" SIZE="ANY">
              <ALTERNATE KIND="CUSTOM" />
              %%Keyword: FIXOrb <advanced>
              <HELP>
              Freezes a subset of VB orbitals (i1, i2, ...).
              </HELP>
              </KEYWORD>

:kword:`FIXStruc`
  .. index::
     single: FIXSTRUC

  ::

    FIXSTRUC
    i1 i2 ...

  Freezes the coefficients for structures :math:`i_1, i_2 \ldots`. Alternatively
  the special keywords :kword:`ALL` or :kword:`NONE` may be used. The
  structures are eliminated from the optimization procedure, but may
  still be affected by normalization or any symmetry keywords present.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="FIXSTRUC" APPEAR="Freeze coefficients" LEVEL="ADVANCED" KIND="STRING">
              %%Keyword: FIXStruc <advanced>
              <HELP>
              Freezes a subset of structure coefficients (i1, i2, ...).
              </HELP>
              </KEYWORD>

:kword:`DELStruc`
  .. index::
     single: DELSTRUC

  ::

    DELSTRUC
    i1 i2 ...

  Deletes the specified structures from the wavefunction. The
  special keywords :kword:`ALL` or :kword:`NONE` may be used. This specification should be compatible
  with the other structure constraints present, as defined by :kword:`SYMELM` and :kword:`ORBREL`.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="DELSTRUC" APPEAR="Delete structures" LEVEL="ADVANCED" KIND="STRING">
              %%Keyword: DELStruc <advanced>
              <HELP>
              Deletes a subset of structures from the wavefunction (i1, i2, ...). Other
              possible values: ALL or NONE.
              </HELP>
              </KEYWORD>

:kword:`ORTHcon`
  .. index::
     single: ORTHCON

  ::

    ORTHCON
    key-1 ...
    key-2 ...
    ...

  The :kword:`ORTHCON` keyword initiates the input of orthogonality
  constraints between pairs/groups of valence bond orbitals.
  The sub-keywords ``key-i`` can be any of :kword:`ORTH`, :kword:`PAIRS`,
  :kword:`GROUP`, :kword:`STRONG` or :kword:`FULL`. Orthogonality constraints
  should be used with discretion. Note that orthogonality constraints
  for an orbital generated from another by symmetry operations (using the
  :kword:`ORBREL` keyword) cannot in general be satisfied. The :kword:`ENDORTHcon`
  keyword can be used to terminate the input of orthogonality constraints. ::

    ORTH
    i1 i2 ...

  Specifies a list of orbitals to be orthogonalized. All overlaps
  between pairs of orbitals in the list are set to zero. ::

    PAIRS i1 i2 ...

  Specifies a simple list of orthogonalization pairs. Orbital :math:`i_1` is
  made orthogonal to :math:`i_2`, :math:`i_3` to :math:`i_4`, etc. ::

    GROUP label i1 i2 ...

  Defines an orbital group to be used with the :kword:`ORTH` or
  :kword:`PAIRS` keyword. The group is referred to by ``label`` which
  can be any three characters beginning with a letter a--z. Labels
  defining different groups can be used together or in combination
  with orbital numbers in :kword:`ORTH` or :kword:`PAIRS`.
  :math:`i_1, i_2 \ldots` specifies
  the list of orbitals in the group. Thus the combination
  :kword:`GROUP` AAA 1 2 :kword:`GROUP` BBB 3 4 :kword:`ORTH` AAA BBB will orthogonalize
  the pairs of orbitals 1--3, 1--4, 2--3 and 2--4. ::

    STRONG

  This keyword is short-hand for strong orthogonality. The only allowed
  non-zero overlaps are between pairs of orbitals (:math:`2n-1`, :math:`2n`). ::

    FULL

  This keyword is short-hand for full orthogonality and is mainly
  useful for testing purposes.

  .. xmldoc:: <GROUP MODULE="CASVB" NAME="ORTHCON" APPEAR="Orthogonality constraints" KIND="BLOCK" LEVEL="ADVANCED">

  .. xmldoc:: %%Keyword: ORTHcon <advanced>
              <HELP>
              Initiates input of orthogonality constraints information.
              Sub-keywords are ORTH, PAIRS, GROUP, STRONG and FULL, as described in the
              manual. The ENDORTHcon keyword terminates the ORTHcon input.
              </HELP>

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="ORTH" APPEAR="Orbitals" LEVEL="ADVANCED" KIND="UNKNOWN" />

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="PAIRS" APPEAR="Pairs" LEVEL="ADVANCED" KIND="UNKNOWN" />

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="GROUP" APPEAR="Group" LEVEL="ADVANCED" KIND="UNKNOWN" />

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="STRONG" APPEAR="Strong orthogonality" LEVEL="ADVANCED" KIND="SINGLE" EXCLUSIVE="FULL" />

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="FULL" APPEAR="Full orthogonality" LEVEL="ADVANCED" KIND="SINGLE" EXCLUSIVE="STRONG" />

  .. xmldoc:: </GROUP>

Optional keywords for wavefunction analysis

.. class:: keywordlist

:kword:`CIWEights`
  .. index::
     single: CIWEIGHTS

  For further details regarding the calculation of weights in :program:`CASVB`, see
  :cite:`casvb7`. ::

    CIWEIGHTS
    key-1 key-2 ... [Nconf]

  Prints weights of the CASSCF wavefunction transformed
  to the basis of nonorthogonal VB structures. For the ``key-i`` options
  see :kword:`VBWEIGHTS` below. Note that the evaluation of inverse overlap
  weights involves an extensive computational overhead for large active
  spaces. Weights are given for the
  total CASSCF wavefunction, as well as the orthogonal complement to
  :math:`\Psi_{\text{VB}}`. The default for the number of configurations requested,
  :math:`N_{\text{conf}}`, is 10. If :math:`N_{\text{conf}} = -1` all configurations are
  included.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="CIWEIGHTS" APPEAR="Print CI weights" LEVEL="ADVANCED" KIND="STRING">
              %%Keyword: CIWEights <advanced>
              <HELP>
              Prints weights of VB structures in the CASSCF wavefunction. Options are
              the same as for VBWEights.
              </HELP>
              </KEYWORD>

:kword:`REPOrt`
  .. index::
     single: REPORT

  ::

    REPORT
    [...
    ENDREPORT]

  Outputs orbital/structure coefficients and derived information.
  The :kword:`ENDREPOrt` keyword can be used to mark the end of the specification
  of a report step.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="REPORT" APPEAR="Report" LEVEL="ADVANCED" KIND="SINGLE">
              %%Keyword: REPOrt <advanced>
              <HELP>
              Outputs orbital/structure coefficients and derived information. Terminated
              by ENDREPOrt.
              </HELP>
              </KEYWORD>

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="ENDREPORT" LEVEL="UNDOCUMENTED" KIND="SINGLE" />

:kword:`(NO)SCORr`
  .. index::
     single: SCORR
     single: NOSCORR

  ::

    (NO)SCORR

  With this option, expectation values of the spin operators
  :math:`(\hat{s}_\mu+\hat{s}_\nu)^2` are evaluated for all pairs of :math:`\mu` and
  :math:`\nu`. Default is :kword:`NOSCORR`. The procedure is described in
  :cite:`casvb8,casvb9,casvb10`.

  This analysis is currently only implemented for spin-coupled wavefunctions.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="SCORR" APPEAR="Spin values" LEVEL="ADVANCED" KIND="SINGLE" EXCLUSIVE="NOSCORR">
              %%Keyword: SCORr <advanced>
              <HELP>
              Performs spin-correlation analysis. Only implemented for spin-coupled
              wavefunctions
              </HELP>
              </KEYWORD>

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="NOSCORR" APPEAR="No spin values" LEVEL="ADVANCED" KIND="SINGLE" EXCLUSIVE="SCORR">
              %%Keyword: NOSCORr <advanced>
              <HELP>
              Disables SCORr
              </HELP>
              </KEYWORD>

:kword:`VBWEights`
  .. index::
     single: VBWEIGHTS

  For further details regarding the calculation of weights in :program:`CASVB`, see
  :cite:`casvb7`. ::

    VBWEIGHTS
    key-1 key-2 ...

  Calculates and outputs weights of the structures in the valence bond
  wavefunction :math:`\Psi_{\text{VB}}`. ``key-i`` specifies the definition of
  nonorthogonal weights to be used, and can be one of:

  :kword:`CHIRGWIN`
    Evaluates Chirgwin--Coulson weights (see :cite:`casvb11`).
  :kword:`LOWDIN`
    Performs a symmetric orthogonalization of the
    structures and outputs the subsequent weights.
  :kword:`INVERSE`
    Outputs "inverse overlap populations" as in :cite:`casvb12`.
  :kword:`ALL`
    All of the above.
  :kword:`NONE`
    Suspends calculation of structure weights.

  The commands :kword:`LOWDIN` and :kword:`INVERSE` require the overlap matrix
  between valence bond structures, so that some additional computational
  overhead is involved.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="VBWEIGHTS" APPEAR="Print VB weights" LEVEL="ADVANCED" KIND="STRING">
              %%Keyword: VBWEights <advanced>
              <HELP>
              Prints weights of VB structures.
              Possible values CHIRGWIN, LOWDIN, INVERSE, ALL and NONE.
              </HELP>
              </KEYWORD>

Optional keywords for further general options

.. class:: keywordlist

:kword:`PREC`
  .. index::
     single: PREC

  ::

    PREC
    iprec iwidth

  Adjusts the precision for printed quantities. In most cases, ``iprec`` simply refers
  to the number of significant digits after the decimal point. Default is ``iprec``\=+8.
  ``iwidth`` specifics the maximum width of printed output, used when determining
  the format for printing arrays.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="PREC" APPEAR="Print precision" LEVEL="ADVANCED" KIND="INTS" SIZE="2" DEFAULT_VALUES="8,110" MIN_VALUE="0">
              %%Keyword: PREC <basic>
              <HELP>
              Adjusts the precision for printed quantities. See manual.
              </HELP>
              </KEYWORD>

:kword:`PRINt`
  .. index::
     single: PRINT

  ::

    PRINT
    i1 i2 ...

  Each number specifies the level of output required at various stages of the execution, according to the
  following convention:

  .. container:: list

    **-1** No output except serious, or fatal, error messages.

    **0**  Minimal output.

    **1**  Standard level of output.

    **2**  Extra output.

  The areas for which output can be controlled are:
  :math:`i_1`

  .. container:: list

    :math:`i_1` Print of input parameters, wavefunction definitions, etc.

    :math:`i_2` Print of information associated with symmetry constraints.

    :math:`i_3` General convergence progress.

    :math:`i_4` Progress of the 2nd-order optimization procedure.

    :math:`i_5` Print of converged solution and analysis.

    :math:`i_6` Progress of variational optimization.

    :math:`i_7` File usage.

  For all, the default output level is +1. If :math:`i_5 \geq 2` VB orbitals will
  be printed in the AO basis (provided that the definition of MOs is
  available).

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="PRINTLEVELS" APPEAR="Print levels" LEVEL="BASIC" KIND="INTS" SIZE="7" DEFAULT_VALUES="1,1,1,1,1,1,1" MIN_VALUE="-1" MAX_VALUE="2">
              %%Keyword: PRINT <basic>
              <HELP>
              Controls the amount of output. See manual.
              </HELP>
              </KEYWORD>

:kword:`SHSTruc`
  Prints overlap and Hamiltonian matrices between VB structures.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="SHSTRUC" APPEAR="Print matrices" LEVEL="ADVANCED" KIND="SINGLE">
              %%Keyword: SHSTruc <advanced>
              <HELP>
              Prints overlap and Hamiltonian matrices between VB structures.
              </HELP>
              </KEYWORD>

:kword:`STATs`
  .. index::
     single: STATS

  ::

    STATS

  Prints timing and usage statistics.

  .. xmldoc:: <KEYWORD MODULE="CASVB" NAME="STATS" APPEAR="Print statistics" LEVEL="BASIC" KIND="SINGLE">
              %%Keyword: STATs <basic>
              <HELP>
              Prints timing and usage statistics.
              </HELP>
              </KEYWORD>

Input example
.............

.. extractfile:: ug/CASVB.input

  &seward
  symmetry
  x y
  basis set
  c.sto-3g....
  c 0 0 -0.190085345
  end of basis
  basis set
  h.sto-3g....
  h 0 1.645045225 1.132564974
  end of basis
  &scf
  occupied
  3 0 1 0
  &rasscf
  inactive
  1 0 0 0
  ras2
  3 1 2 0
  nactel
  6 0 0
  lumorb
  &casvb

.. index::
   single: CASVB; Plotting

Viewing and plotting VB orbitals
................................

In many cases it can be helpful to view the shape of the converged valence bond orbitals, and
Molcas therefore provides two facilities for doing this. For the Molden program, an interface file
is generated at the end of each :program:`CASVB` run (see also :numref:`UG:sec:Molden`).
Alternatively a :program:`CASVB` run may be followed by :program:`RASSCF` to get orbitals
(:numref:`UG:sec:rasscf`) and :program:`GRID_IT` with the :kword:`VB` specification
(:numref:`UG:sec:gridit`), in order to generate a three-dimensional grid, for viewing, for example,
with :program:`LUSCUS` program.

.. xmldoc:: </MODULE>
