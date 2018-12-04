.. index::
   single: Program; CHCC
   single: CHCC

.. _sec\:chcc:

:program:`chcc`
===============

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="CHCC">
            %%Description:
            <HELP>
            The CHCC is a Closed-Shell Coupled-Clusters Singles and Doubles program
            based exclusively on the Cholesky (or RI) decomposed 2-electron integrals aimed
            towards calculation of large system. Use of point-group symmetry is not
            implemented.
            It requires RUNFILE from the SCF module and Cholesky vectors from SEWARD.
            </HELP>

The :program:`CHCC` is a Closed-Shell Coupled-Clusters Singles and Doubles
program based exclusively on the Cholesky (or RI) decomposed 2-electron integrals
aimed towards calculation of large systems on highly parallel architectures. Use of
point-group symmetry is not implemented. Main advantage compared to the
:program:`CCSDT` module in |molcas| is in its more efficient parallelization and
dramatically lowered memory (and eventually disk) requirements.

.. For further details the reader is referred to the section :ref:`TUT:sec:chcc`.

.. index::
   pair: Dependencies; CHCC

.. _sec\:chcc_dependencies:

Dependencies
------------

:program:`CHCC` requires a previous run of the RHF :program:`SCF` program
to produce molecular orbitals and orbital energies stored in :file:`RUNFILE`.
The :program:`SCF` program (as well as :program:`SEWARD`) must be run
in Cholesky/RI mode.

The algorithm used for almost complete elimination of the :program:`CHCC`
limits in calculated system size due to the computer memory bottleneck relies
on blocking of the virtual orbitals. Number of blocks (further also referred to as the
"large" segmentation, :kword:`LARGe`), :math:`N'`, should be as small as
possible, because increasing of the segmentation brings in more CPU and I/O overhead.
Furthermore, blocking can be "fine tuned" by, so called, "small" segmentation (:kword:`SMALl`), :math:`N''`,
which affects only the (typically) most demanding :math:`O^2V^4` scaling
terms. The "large" segmentation can range from 1 to 32, "small" segmentation from 1 to 8, but
their product, i.e. "large" |x| "small" must be no more than 64.

Selected blocking also determines the
number of "independent" parallel tasks that must be executed in each iteration of
the CCSD equations. In other words, particular segmentation predetermines the optimal
number of computational nodes (i.e., if the best possible parallelization is desired).
If the requested "large" segmentation is :math:`N'`, then :math:`N'^2` terms scaling as
:math:`O^3V^3` and :math:`N'^2/2` terms scaling as :math:`O^2V^4` result.
Depending on which of these terms dominated in the calculations (:math:`O^3V^3`
is more demanding for systems with large number of occupied orbitals and rather small
basis set, while :math:`O^2V^4` dominated for relatively large basis sets,
i.e. large number of virtual orbitals), number of these task should be divisible by the number of
computational nodes for optimal performance. To make it simple, as a rule of thumb, :math:`N'^2/2`
should be divisible by the number of nodes, since the :math:`O^3V^3` are typically twice less
expensive then the :math:`O^2V^4` step. Otherwise, any reasonable (i.e. the number
of tasks is larger than the number of computational nodes, obviously) combination is allowed.

.. index::
   pair: Files; CHCC

.. _sec\:chcc_files:

Files
-----

Input files
...........

:program:`CHCC` will use the following input
files: :file:`CHVEC`, :file:`CHRED`, :file:`CHORST`, :file:`RUNFILE`,
and :file:`CHOR2F`
(for more information see :ref:`UG:sec:files_list`).

Output files
............

.. class:: filelist

:file:`L0xxxx`, :file:`L1xxxx`, :file:`L2xxxx`
  MO-transformed Cholesky vectors

:file:`T2xxxx`
  T2 :math:`(ij,a'b')` excitation amplitudes

:file:`RstFil`
  Communication file containing T1 amplitudes, restart informations, etc.

.. index::
   pair: Input; CHCC

.. _sec\:chcc_input:

Input
-----

The input for each module is preceded by its name like: ::

  &CHCC

Optional keywords

.. class:: keywordlist

:kword:`TITLe`
  This keyword is followed by one title line.

  .. xmldoc:: <KEYWORD MODULE="CHCC" NAME="TITLE" APPEAR="Title" KIND="STRINGS" SIZE="10" LEVEL="BASIC">
              %%Keyword: TITLe <basic>
              <HELP>
              Enter up to ten title lines. Do not put any keyword in the beginning of a title line.
              </HELP>
              </KEYWORD>

:kword:`FROZen`
  Integer on the following line specifies number of inactive occupied
  orbitals in the CCSD calculation. (Default=0)

  .. xmldoc:: <KEYWORD MODULE="CHCC" NAME="FROZ" APPEAR="Frozen orbitals" KIND="INT" LEVEL="BASIC" MIN_VALUE="0" DEFAULT_VALUE="0">
              %%Keyword: FROZen <basic>
              <HELP>
              Specifies number of inactive occupied orbitals in the CCSD procedure
              </HELP>
              </KEYWORD>

:kword:`DELEted`
  Integer on the following line specifies number of inactive virtual
  orbitals in the CCSD calculation. (Default=0)

  .. xmldoc:: <KEYWORD MODULE="CHCC" NAME="DELE" APPEAR="Deleted orbitals" KIND="INT" LEVEL="BASIC" MIN_VALUE="0" DEFAULT_VALUE="0">
              %%Keyword: DELEted <basic>
              <HELP>
              Specifies number of inactive virtual orbitals in the CCSD procedure
              </HELP>
              </KEYWORD>

:kword:`LARGe`
  Integer on the following line specifies the main segmentation of the virtual orbitals.
  Value must be between 1 (no segmentation) and 32. Product of Large and Small segmentation
  must be lower than 64. (Default=1)

  .. xmldoc:: <KEYWORD MODULE="CHCC" NAME="LARG" APPEAR="Large segmentation" KIND="INT" LEVEL="BASIC" MIN_VALUE="1" MAX_VALUE="32" DEFAULT_VALUE="1">
              %%Keyword: LARGe <basic>
              <HELP>
              Specifies the segmentation of virtual orbitals
              </HELP>
              </KEYWORD>

:kword:`SMALl`
  Integer on the following line specifies the auxiliary segmentation of the virtual orbitals.
  Value must be between 1 (no segmentation) and 8. Product of Large and Small segmentation
  must be lower than 64. Small segmentation doesn't generate extra parallel tasks.
  (Default=1)

  .. xmldoc:: <KEYWORD MODULE="CHCC" NAME="SMAL" APPEAR="Small segmentation" KIND="INT" LEVEL="BASIC" MIN_VALUE="1" MAX_VALUE="8" DEFAULT_VALUE="1">
              %%Keyword: SMALl <basic>
              <HELP>
              Specifies the auxiliary segmentation of virtual orbitals
              </HELP>
              </KEYWORD>

:kword:`CHSEgmentation`
  Integer on the following line specifies the block size of the auxiliary (Cholesky/RI)
  index. Value must be lower than the minimal dimension of the auxiliary index on each
  computational node. (Default=100)

  .. xmldoc:: <KEYWORD MODULE="CHCC" NAME="CHSE" APPEAR="Auxiliary block size" KIND="INT" LEVEL="BASIC" MIN_VALUE="1" DEFAULT_VALUE="100">
              %%Keyword: CHSEgmentation <basic>
              <HELP>
              Specifies the block size of auxiliary (Cholesky/RI) index
              </HELP>
              </KEYWORD>

:kword:`MHKEy`
  Integer on the following line specifies if library BLAS (MHKEy=1) or hard-coded
  fortran vector-vector, matrix-vector and matrix-matrix manipulation is used.
  (Default=1)

  .. xmldoc:: <KEYWORD MODULE="CHCC" NAME="MHKE" APPEAR="Use BLAS" KIND="INT" LEVEL="BASIC" MIN_VALUE="0" MAX_VALUE="1" DEFAULT_VALUE="1">
              %%Keyword: MHKEy <basic>
              <HELP>
              Specifies if BLAS libraries (=1) or hard-code fortran is used.
              </HELP>
              </KEYWORD>

:kword:`NOGEnerate`
  This keyword specifies that the pre-CCSD steps (regeneration of integrals from
  the Cholesky/RI vectors, etc.) are skipped.
  (Default=OFF)

  .. xmldoc:: <KEYWORD MODULE="CHCC" NAME="NOGE" APPEAR="Skip pre-CCSD steps" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: NOGEnerate <basic>
              <HELP>
              Pre-CCSD steps, like integrals generation, etc. are skipped.
              </HELP>
              </KEYWORD>

:kword:`ONTHefly`
  This keyword specifies that all integral types scaling steeper then :math:`O^2V^2`
  are generated "on-the-fly" from the Cholesky/RI vectors. Use of this keyword leads
  to dramatically savings of the disk resources, but leads to significant arithmetic
  overhead. Keywords "ONTHefly"
  and "PRECalculate" are mutually exclusive.
  (Default=OFF)

  .. xmldoc:: <SELECT MODULE="CHCC" NAME="INTEGRALS" APPEAR="Integrals" CONTAINS="ONTH,PREC">

  .. xmldoc:: <KEYWORD MODULE="CHCC" NAME="ONTH" APPEAR="On the fly" KIND="SINGLE" LEVEL="BASIC" EXCLUSIVE="PREC">
              %%Keyword: ONTHefly <basic>
              <HELP>
              Integrals with 3- and 4-virtual indexes are generated "on-the-fly".
              </HELP>
              </KEYWORD>

:kword:`PRECalculate`
  This keyword specifies that all integral are precalculated before the
  CCSD iterative procedure starts. Use of this keyword leads to significant
  consumption of the disk space, especially is single-processor runs.
  (Default=ON)

  .. xmldoc:: <KEYWORD MODULE="CHCC" NAME="PREC" APPEAR="Precalculate (default)" KIND="SINGLE" LEVEL="BASIC" EXCLUSIVE="ONTH">
              %%Keyword: PRECalculate <basic>
              <HELP>
              All integrals are precalculated prior to the CCSD iterations.
              </HELP>
              </KEYWORD>

  .. xmldoc:: </SELECT>

:kword:`NODIstribute`
  This keyword (in combination with the "PRECalculate" keyword) specifies that all
  integral are stored on each computational node. In case of all integrals being
  stored on each node, extra permutation symmetry can be applied, thus leading to
  significant savings of the disk space. However, in case of massively parallel runs
  (i.e. more than ~8 nodes), savings from keeping only subset of integrals
  required on particular node are more significant than savings due to permutational
  symmetry. (Default=OFF)

  .. xmldoc:: <KEYWORD MODULE="CHCC" NAME="NODI" APPEAR="No distribute" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: NODIstribute <basic>
              <HELP>
              All integrals are precalculated on all computational nodes.
              </HELP>
              </KEYWORD>

:kword:`JOINlkey`
  The parameter on the following line specifies, which algorithm is used for
  precalculation and of the integrals in parallel run. In parallel runs, :program:`SEWARD`
  produces AO Cholesky/RI vectors segmented in auxiliary index over
  parallel nodes. Depending on the network bandwidth and computational power
  of each node, different algorithms can lead to optimal performance.
  Following options are available:

  .. container:: list

    0 --- None: no cumulation of Cholesky/RI vectors is needed (debug only).

    1 --- Minimal: Cholesky/RI vectors are cumulated prior to integral precalculation. Low network bandwidth is required.

    2 --- Medium: :math:`O^2V^2` integrals are generated from local Cholesky/RI vectors and cumulated along with the Cholesky/RI vectors afterwards.
    Other integrals are calculated from cumulated intermediates.

    3 --- Full: All integrals are generated from local Cholesky/RI vectors and cumulated afterwards. High network bandwidth is required.

  (Default=2)

  .. xmldoc:: <KEYWORD MODULE="CHCC" NAME="JOIN" APPEAR="Parallel integral generation" KIND="CHOICE" LIST="0: None (debug),1: Minimal,2: Medium,3: Full" LEVEL="BASIC" DEFAULT_VALUE="2">
              <HELP>
              Choose the type of parallel integral generation
              </HELP></KEYWORD>
              %%Keyword: JOINlkey <basic>
              Sets the type of parallel integral generation
              ||0 - None (debug)
              ||1 - Minimal (low network bandwdith required)
              ||2 - Medium
              ||3 - Full (high network bandwidth required)

:kword:`MAXIterations`
  Integer on the following line specifies maximum number of CCSD iteration
  (Default=40)

  .. xmldoc:: <KEYWORD MODULE="CHCC" NAME="MAXI" APPEAR="Maximum iterations" KIND="INT" LEVEL="BASIC" MIN_VALUE="0" DEFAULT_VALUE="40">
              %%Keyword: MAXIterations <basic>
              <HELP>
              Maximum number of CCSD iterations.
              </HELP>
              </KEYWORD>

:kword:`RESTart`
  This keyword specifies that CCSD calculation is restarted from previous
  run. This keyword is currently under development,
  thus disabled. (Default=OFF)

  .. xmldoc:: <KEYWORD MODULE="CHCC" NAME="REST" APPEAR="Restart" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: RESTart <basic>
              <HELP>
              Restart from previous run. Currently disabled.
              </HELP>
              </KEYWORD>

:kword:`THREshold`
  Double precision floating point number on the following line specifies
  the convergence threshold for the CCSD correlation energy.
  (Default=1.0d-6)

  .. xmldoc:: <KEYWORD MODULE="CHCC" NAME="THRE" APPEAR="Convergence threshold" KIND="REAL" LEVEL="BASIC" MIN_VALUE="0.0" DEFAULT_VALUE="1.0d-6">
              %%Keyword: THREshold <basic>
              <HELP>
              Convergence threshold for the CCSD correlation energy.
              </HELP>
              </KEYWORD>

:kword:`PRINtkey`
  The integer on the following line specifies the print level in output

  .. container:: list

    1  --- Minimal

    2  --- Minimal + timings of each step of the CCSD iterations

    10 --- Debug

  (Default=1)

  .. xmldoc:: <KEYWORD MODULE="CHCC" NAME="PRIN" APPEAR="Print level" KIND="CHOICE" LIST="1: Minimal,2: Minimal + timings,10: Debug" LEVEL="ADVANCED" DEFAULT_VALUE="1">
              <HELP>
              Choose the print level
              </HELP>
              </KEYWORD>
              %%Keyword: PRINtkey <advanced>
              Sets the print level
              || 1 - Minimal
              || 2 - Minimal + timings
              ||10 - Debug

:kword:`END of input`
  This keyword indicates that there is no more input
  to be read.

::

  &CHCC &END
  Title
  Benzene dimer
  Frozen
  12
  Deleted
  0
  Large
  4
  Small
  2
  CHSEgment
  100
  Precalculate
  Join
  2
  Maxiter
  50
  Threshold
  1.0d-6
  Print
  2
  End of Input

.. xmldoc:: </MODULE>
