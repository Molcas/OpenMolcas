.. index::
   single: Program; Surfacehop
   single: Surfacehop

.. _UG\:sec\:surfacehop:

:program:`surfacehop`
=====================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="SURFACEHOP">
            %%Description:
            <HELP>

            QUICK DESCRIPTION OF THE PROGRAM

            </HELP>

This module deals with surface hop semiclassical molecular dynamics (SHMD) and has to be used together with module DYNAMIX. Its purpose is the calculation of the relax root for the next step of the SHMD. In this moment the implemented algorithm under this module is the Tully's fewest switches :cite:`Tully1990`, along with the Hammes-Schiffer/Tully scheme :cite:`Hammes-Schiffer1994` and the decoherence correction proposed by Granucci and Persico :cite:`Granucci2007`.

.. _UG\:sec\:surfacehop_output_files:

Output files
............

.. class:: filelist

:file:`RUNFILE`
  Surface hop information such as Amatrix and CI coefficients for previous steps are stored in this file.

.. _UG\:sec\:surfacehop_inp:

Input
-----

::

  &Gateway
  coord=$Project.xyz
  basis=6-31G*
  group=nosym

  >> EXPORT MOLCAS_MAXITER=400
  >> DOWHILE

  &Seward

  &rasscf
   jobiph
   cirestart
   nactel = 6 0 0
   inactive = 23
   ras2 = 6
   ciroot = 2 2 1
   prwf = 0.0
   mdrlxroot = 2

  &Surfacehop
   tully
   decoherence = 0.1
   psub

  &alaska

  &Dynamix
   velver
   dt = 41.3
   velo = 1
   thermo = 0
  >>> End Do

General keywords
................

.. class:: keywordlist

:kword:`TULLY`
  This keyword enables the Tully--Hammes-Schiffer integration of the TDSE for the Tully Surface Hop Algorithm. If you use this keyword you should not use the :kword:`HOP` keyword in :program:`DYNAMIX`.

  .. xmldoc:: <KEYWORD MODULE="SURFACEHOP" NAME="TULLY" APPEAR="Tully surface hop" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: TULLy <basic>
              <HELP>
              This keyword enables the Tully--Hammes-Schiffer integration of the TDSE for the Tully Surface Hop Algorithm.
              </HELP>
              </KEYWORD>

:kword:`DECOHERENCE`
  This keyword must be used after the :kword:`TULLY` keyword. It enables the decoherence correction in the population density matrix as reported by Persico--Granucci. The value is called decay factor and it is usually 0.1 hartree. It can be seen as how strongly this correction is applied. It is recommendable to leave it to 0.1, unless you really know what your're doing.

  .. xmldoc:: <KEYWORD MODULE="SURFACEHOP" NAME="DECOHERENCE" APPEAR="Decoherence correction" KIND="REAL" LEVEL="ADVANCED" DEFAULT_VALUE="0.1" REQUIRE="TULLY">
              %%Keyword: DECOherence <advanced>
              This keyword must be used after the TULLY keyword.
              <HELP>
              It enables the decoherence correction in the population density matrix as reported by Persico--Granucci.
              </HELP>
              </KEYWORD>

:kword:`SUBSTEP`
  This keyword must be used after the :kword:`TULLY` keyword. This keyword specifies how many steps of integration we use to interpolate/extrapolate between two Newton's consecutive steps. The default is usually a good compromise between quickness and precision (200 substeps each femtoseconds of MD).

  .. xmldoc:: <KEYWORD MODULE="SURFACEHOP" NAME="SUBSTEP" APPEAR="Electronic integration substeps" KIND="INT" LEVEL="ADVANCED" DEFAULT_VALUE="200" MIN_VALUE="0" REQUIRE="TULLY">
              %%Keyword: SUBStep <advanced>
              This keyword must be used after the TULLY keyword.
              <HELP>
              It specifies how many steps of integration we use to interpolate/extrapolate between two Newton's consecutive steps.
              </HELP>
              </KEYWORD>

:kword:`PSUB`
  This keyword must be used after the :kword:`TULLY` keyword. To print in |molcas| output :math:`\mat{D}` matrix, :math:`\mat{A}` matrix, :math:`\mat{B}` matrix, probabilities, randoms, population and energies at each substep (quite verbose, but gives you a lot of useful information).

  .. xmldoc:: <KEYWORD MODULE="SURFACEHOP" NAME="PSUB" APPEAR="Verbose output for each substep" KIND="SINGLE" LEVEL="BASIC" REQUIRE="TULLY">
              %%Keyword: PSUB <basic>
              This keyword must be used after the TULLY keyword.
              <HELP>
              To print in molcas output D matrix, A matrix, B matrix, probabilities, randoms, population and energies at each substep.
              </HELP>
              </KEYWORD>

:kword:`DMTX`
  This keyword must be used after the :kword:`TULLY` keyword. With this keyword you can start your calculation with an initial :math:`\mat{A}` matrix (population density matrix). It is a complex matrix. In the first line after the keyword you must specify its dimension :math:`N`. Then :math:`N` lines (:math:`N` values each line) with the real part of the matrix followed by :math:`N` more lines with the imaginary part.

  .. xmldoc:: <KEYWORD MODULE="SURFACEHOP" NAME="DMTX" APPEAR="Initial population density matrix" KIND="UNKNOWN" LEVEL="ADVANCED" REQUIRE="TULLY">
              %%Keyword: DMTX <advanced>
              This keyword must be used after the TULLY keyword.
              <HELP>
              Initial A matrix (population density matrix). It is a complex matrix.
              In the first line after the keyword you must specify its dimension N. Then N lines (with N values each line) with the REAL part of the matrix followed by N more lines with the imaginary part.
              </HELP>
              </KEYWORD>

:kword:`FRANDOM`
  This keyword must be used after the :kword:`TULLY` keyword. It fixes the random number to one provided by the user, in case a deterministic trajectory is needed

  .. xmldoc:: <KEYWORD MODULE="SURFACEHOP" NAME="FRANDOM" APPEAR="Random number constant (deterministic MD)" KIND="REAL" LEVEL="ADVANCED" REQUIRE="TULLY">
              %%Keyword: FRANdom <advanced>
              This keyword must be used after the TULLY keyword.
              <HELP>
              It fixes the random number to one provided by the user, in case a deterministic trajectory is needed.
              </HELP>
              </KEYWORD>

:kword:`ISEED`
  This keyword must be used after the :kword:`TULLY` keyword. The initial seed number is read from the input file. Then, seed numbers are modified (in a deterministic way), saved in the :file:`RunFile` and read in the next call to the module. This way, MD simulations are reproducible.

  .. xmldoc:: <KEYWORD MODULE="SURFACEHOP" NAME="ISEED" APPEAR="Initial seed number (reproducible MD)" KIND="INT" LEVEL="ADVANCED" REQUIRE="TULLY">
              %%Keyword: ISEEd <advanced>
              This keyword must be used after the TULLY keyword.
              <HELP>
              The initial seed number is read from the input file.
              Then, seed numbers are modified (in a deterministic way), saved in the RunFile and read in the next call to the module.
              This way, MD simulations are reproducible.
              </HELP>
              </KEYWORD>

:kword:`MAXHOP`
  This keyword must be used after the :kword:`TULLY` keyword. It specifies how many non-adiabatic transitions are allowed between electronic states.

  .. xmldoc:: <KEYWORD MODULE="SURFACEHOP" NAME="MAXHOP" APPEAR="Maximum number of hops allowed" KIND="INT" LEVEL="ADVANCED" REQUIRE="TULLY">
              %%Keyword: MAXHop <advanced>
              This keyword must be used after the TULLY keyword.
              <HELP>
              It specifies how many non-adiabatic transitions are allowed between electronic states.
              </HELP>
              </KEYWORD>

:kword:`H5RESTART`
  This keyword allows to restart a surface hopping trajectory calculation from an HDF5 file.
  The name of the restart file is given on the next line.

  .. xmldoc:: <KEYWORD MODULE="SURFACEHOP" NAME="H5RESTART" APPEAR="Restart the surface hopping trajectory from an H5 file" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: H5REstart <advanced>
              <HELP>
              Restarts a surface hopping trajectory calculation from an HDF5 file, whose name is given on the next line.
              </HELP>
              </KEYWORD>

Input examples
..............

This example shows an excited state CASSCF MD simulation
of a methaniminium cation using the Tully Surface Hop algorithm.
Within the :program:`Surfacehop` module The keyword :kword:`TULLY` enables the TDSE integration. The options used in this case are:
(:kword:`SUBSTEP`\=200) to specify 200 substep of electronic integration between Newton's,
(:kword:`DECOHERENCE`\=1) to deal with the decoherence using a decay constant of 0.1 hartree and
(:kword:`PSUB`) to print the substeps matrices verbosely into the |molcas| log.

.. extractfile:: ug/surfacehopTULLY.input

  &GATEWAY
   COORD
   6
   Angstrom
   C  0.00031448  0.00000000  0.04334060
   N  0.00062994  0.00000000  1.32317716
   H  0.92882820  0.00000000 -0.49115611
   H -0.92846597  0.00000000 -0.49069213
   H -0.85725321  0.00000000  1.86103989
   H  0.85877656  0.00000000  1.86062860
   BASIS= 3-21G
   GROUP= nosym

  >> EXPORT MOLCAS_MAXITER=1000
  >> DOWHILE

  &SEWARD

  >> IF ( ITER = 1 )

  &RASSCF
    LUMORB
   FileOrb= $Project.GssOrb
   Symmetry= 1
   Spin= 1
   nActEl= 2 0 0
   Inactive= 7
   RAS2= 2
   CIroot= 3 3 1

  >> COPY $Project.JobIph $Project.JobOld

  >> ENDIF

  &RASSCF
   JOBIPH; CIRESTART
   Symmetry= 1
   Spin= 1
   nActEl= 2 0 0
   Inactive= 7
   RAS2= 2
   CIroot= 3 3 1
   MDRLXR= 2

  >> COPY $Project.JobIph $Project.JobOld

  &surfacehop
   TULLY
   SUBSTEP = 200
   DECOHERENCE = 0.1
   PSUB

  &ALASKA

  &Dynamix
   VELVer
   DT= 10.0
   VELO= 3
   THER= 2
   TEMP=300

  >> END DO

.. xmldoc:: </MODULE>
