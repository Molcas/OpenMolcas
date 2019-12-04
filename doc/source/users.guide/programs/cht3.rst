.. index::
   single: Program; CHT3
   single: CHT3

.. _sec\:cht3:

:program:`cht3`
===============

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="CHT3">
            %%Description:
            <HELP>
            CHT3 is a Closed-Shell Coupled-Clusters perturbative triples
            program based exclusively on the Cholesky (or RI) decomposed 2-electron integrals
            aimed towards calculation of large systems on highly parallel architectures. Use of
            point-group symmetry is not implemented.
            It requires RUNFILE from the SCF module, T1 and T2 excitation amplitudes and MO
            transformed Cholesky/RI vectors from CHCC.
            </HELP>

:program:`CHT3` is a Closed-Shell Coupled-Clusters perturbative triples
program based exclusively on the Cholesky (or RI) decomposed 2-electron integrals
aimed towards calculation of large systems on highly parallel architectures. Use of
point-group symmetry is not implemented. Main advantage compared to the
:program:`CCSDT` module is |molcas| is in its more efficient parallelization and
dramatically lowered memory (and eventually disk) requirements.

.. For further details the reader is referred to :numref:`TUT:sec:cht3`.

.. index::
   pair: Dependencies; CHT3

.. _sec\:cht3_dependencies:

Dependencies
------------

:program:`CHT3` requires previous run of the :program:`CHCC` Cholesky/RI
based CCSD program to produce T1 and T2 excitation amplitudes stored in :file:`T2xxxx`
and :file:`RstFil` files.
The :program:`CHCC` program (as well as :program:`SEWARD` and :program:`SCF`) must be run
in Cholesky/RI mode.

The algorithm used for almost complete elimination of the :program:`CHT3`
limits in calculated system size due to the computer memory bottleneck relies
on blocking of the virtual orbitals. Size of blocks is, unlike in :program:`CHCC` program, determined automatically for optimal performance.

.. index::
   pair: Files; CHT3

.. _sec\:cht3_files:

Files
-----

Input files
...........

.. class:: filelist

:file:`RUNFILE`
  File for communication of auxiliary information.

:file:`L0xxxx`, :file:`L1xxxx`, :file:`L2xxxx`
  MO-transformed Cholesky vectors

:file:`T2xxxx`
  T2 :math:`(ij,a'b')` excitation amplitudes

:file:`RstFil`
  Communication file containing T1 amplitudes, restart informations, etc.

Intermediate files
..................

All the intermediate files are created, used and removed
automatically, unless you yourself create a link or a file
with the specified name.

.. class:: filelist

:file:`KMATAA`, :file:`KMATBA`, :file:`LMATAA`, :file:`LMATBA`
  Temporary integral files

Output files
............

None

.. index::
   pair: Input; CHT3

.. _sec\:cht3_input:

Input
-----

The input for each module is preceded by its name like: ::

  &CHT3

.. class:: keywordlist

:kword:`TITLe`
  This keyword starts the reading of title lines,
  with the number of title lines limited to 10.
  Reading the input as title lines is stopped as soon
  as the input parser detects one of the other keywords,
  however only ten lines will be accepted.
  This keyword is *optional*.

  .. xmldoc:: <KEYWORD MODULE="CHT3" NAME="TITLE" APPEAR="Title" KIND="CUSTOM" LEVEL="BASIC">
              %%Keyword: TITLe <basic>
              <HELP>
              Enter up to ten title lines. Do not put any keyword in the beginning of a title line.
              </HELP>
              </KEYWORD>

:kword:`FROZen`
  Integer on the following line specifies number of inactive occupied
  orbitals in the (T) calculation. This keyword is *optional*. (Default=0)

  .. xmldoc:: <KEYWORD MODULE="CHT3" NAME="FROZ" APPEAR="Frozen orbitals" KIND="INT" LEVEL="BASIC" MIN_VALUE="0" DEFAULT_VALUE="0">
              %%Keyword: FROZen <basic>
              <HELP>
              Specifies number of inactive occupied orbitals in the (T) procedure
              </HELP>
              </KEYWORD>

:kword:`DELEted`
  Integer on the following line specifies number of inactive virtual
  orbitals in the (T) calculation. This keyword is *optional*. (Default=0)

  .. xmldoc:: <KEYWORD MODULE="CHT3" NAME="DELE" APPEAR="Deleted orbitals" KIND="INT" LEVEL="BASIC" MIN_VALUE="0" DEFAULT_VALUE="0">
              %%Keyword: DELEted <basic>
              <HELP>
              Specifies number of inactive virtual orbitals in the (T) procedure
              </HELP>
              </KEYWORD>

:kword:`LARGe`
  Integer on the following line specifies the main segmentation of the virtual orbitals
  used in previous CCSD run.

  .. xmldoc:: <KEYWORD MODULE="CHT3" NAME="LARG" APPEAR="Large segmentation" KIND="INT" LEVEL="BASIC" MIN_VALUE="1" MAX_VALUE="32" DEFAULT_VALUE="1">
              %%Keyword: LARGe <basic>
              <HELP>
              Specifies the segmentation of virtual orbitals
              </HELP>
              </KEYWORD>

:kword:`MHKEy`
  Integer on the following line specifies if library BLAS (MHKEy=1) or hard-coded
  fortran vector-vector, matrix-vector and matrix-matrix manipulation is used.
  This keyword is *optional*. (Default=1)

  .. xmldoc:: <KEYWORD MODULE="CHT3" NAME="MHKE" APPEAR="Use BLAS" KIND="INT" LEVEL="BASIC" MIN_VALUE="0" MAX_VALUE="1" DEFAULT_VALUE="1">
              %%Keyword: MHKEy <basic>
              <HELP>
              Specifies if BLAS libraries (=1) or hard-code fortran is used.
              </HELP>
              </KEYWORD>

:kword:`NOGEnerate`
  This keyword specifies that the pre-(T) steps (generation of integrals from
  the Cholesky/RI vectors, etc.) are skipped. This keyword can be used for
  restarting the (T) calculation if the required integrals were already generated.
  This keyword is *optional*. (Default=OFF)

  .. xmldoc:: <KEYWORD MODULE="CHT3" NAME="NOGE" APPEAR="Skip pre-(T) steps" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: NOGEnerate <basic>
              <HELP>
              Pre-(T) steps, like integrals generation, etc. are skipped.
              </HELP>
              </KEYWORD>

:kword:`NOTRiples`
  This keyword specifies that the post integral preparation steps, i.e.
  the real calculation of (T) contribution will not be done. Job can be restarted
  from this point using the :kword:`NOGEnerate` keyword.
  This keyword is *optional*. (Default=OFF)

  .. xmldoc:: <KEYWORD MODULE="CHT3" NAME="NOTR" APPEAR="No triples" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: NOTRiples <basic>
              <HELP>
              Program stops after generation of temporary integral files
              </HELP>
              </KEYWORD>

:kword:`ALOOp`
  Two integers on the following line specify first and last triplet of
  virtual orbitals blocks to be calculated in the first ("A loop") of
  the two parts of the (T) calculation. Using this keyword enables user
  to split the (T) calculation into separate jobs. Information about
  the total number of triplets in the "A loop" can be found in the
  output of the "preparation" step of the (T) program. Values
  -1, -1 mean, that the whole "A loop" is either executed or skipped,
  depending on the parameters of the BLOOp keyword.
  This keyword is *optional*. (Default=-1,-1)

  .. xmldoc:: <KEYWORD MODULE="CHT3" NAME="ALOO" APPEAR="A loop" KIND="INTS" SIZE="2" LEVEL="BASIC" MIN_VALUE="-1" DEFAULT_VALUE="-1">
              %%Keyword: ALOOp <basic>
              <HELP>
              Specifies the range of triplets of virtual orbitals blocks from
              the fist of two parts of (T) program to be calculated.
              </HELP>
              </KEYWORD>

:kword:`BLOOp`
  Two integers on the following line specify first and last triplet of
  virtual orbital block to be calculated in the second ("B loop") of
  two parts of the (T) calculation. Using this keyword enables user
  to split the (T) calculation into separate jobs. Information about
  the total number of triplets in the "B loop" can be found in the
  output of the "preparation" step of the (T) program. Values
  -1, -1 mean, that the whole "B loop" is either executed or skipped,
  depending on the values of the ALOOp keyword.
  This keyword is *optional*. (Default=-1,-1)

  .. xmldoc:: <KEYWORD MODULE="CHT3" NAME="BLOO" APPEAR="B loop" KIND="INTS" SIZE="2" LEVEL="BASIC" MIN_VALUE="-1" DEFAULT_VALUE="-1">
              %%Keyword: BLOOp <basic>
              <HELP>
              Specifies the range of triplets of virtual orbitals blocks from
              the second of two parts of (T) program to be calculated.
              </HELP>
              </KEYWORD>

:kword:`PRINtkey`
  The integer on the following line specifies the print level in output

  .. container:: list

    1  --- Minimal

    2  --- Minimal + timings of each (T) step

    10 --- Debug

  This keyword is *optional*. (Default=1)

  .. xmldoc:: <KEYWORD MODULE="CHT3" NAME="PRIN" APPEAR="Print level" KIND="CHOICE" LIST="1: Minimal,2: Minimal + timings,10: Debug" LEVEL="ADVANCED" DEFAULT_VALUE="1">
              <HELP>
              Choose the print level
              </HELP>
              </KEYWORD>
              %%Keyword: PRINtkey <advanced>
              Sets the print level
              || 1 - Minimal
              || 2 - Minimal + timings
              ||10 - Debug

::

  &CHT3
  Title  = Benzene dimer
  Frozen = 12
  Large  = 4
  ALOOp  = 20 120
  BLoop  = 1 250
  Print  = 2

.. xmldoc:: </MODULE>
