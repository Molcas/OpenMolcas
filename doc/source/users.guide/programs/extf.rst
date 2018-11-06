.. index::
   single: Program; Extf
   single: Extf

.. _UG\:sec\:extf:

:program:`extf`
=====================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="extf">
            %%Description:
            <HELP>

            QUICK DESCRIPTION OF THE PROGRAM

            </HELP>

This module calculates the contribution of an external force that is acting on the system. It applies the modification directly on the gradient and it needs to be called after the execution of :program:`ALASKA`, in an optimization or molecular dynamics calculation. At present time, just the `LINEAR` keyword is present, that applies a constant linear force between two atoms :cite:`valentini2017optomechanical`. The following input example is a semiclassical molecular dynamics with tully surface hop, where a linear force of 2.9 Hartree/Bohr is applied between atom 1 and atom 2.

.. _UG\:sec\:extf_inp:

Input
-----

::

  &Gateway
  coord=$Project.xyz
  basis=6-31G*
  group=nosym

  >> FOREACH A in (1..400)

  &Seward

  &rasscf
   nactel = 6 0 0
   inactive = 23
   ras2 = 6
   ciroot = 2 2 1
   prwf = 0.0
   mdrlxroot = 2

  &alaska

  &surfacehop
   tully
   decoherence = 0.1
   psub

  &Extf
   LINEAR
   1
   2
   2.9
   0

  &Dynamix
   velver
   dt = 41.3
   velo = 1
   thermo = 0
  >>> End Do

General keywords
................

.. class:: keywordlist

:kword:`LINEAR`
  This keyword enables the linear external force.

  .. xmldoc:: <KEYWORD MODULE="extf" NAME="LINEAR" APPEAR="Linear external force" KIND="STRINGS" LEVEL="ADVANCED">
              %%Keyword: LINEar <advanced>
              This keyword enables the linear external force between two atoms.
              <HELP>
              This keyword works by specifying 4 parameters, each one in its own line after the keyword itself. First parameter (Integer) is the first atom number following the numeration of the geometry. Second parameter (Integer) is the second atom number. Third parameter is the force (Float) in Hartree/Bohr applied along the vector between the two atoms. Fourth parameter is 0 or 1 (Bool), where 0 indicates a repulsive force, and 1 is for an attractive force.
              </HELP>
              </KEYWORD>

Input examples
..............

This example shows an excited state CASSCF MD simulation
of a methaniminium cation using the Tully Surface Hop algorithm. In the simulation, the carbon and the nitrogen are pulled apart with a constant force of 1.5 Hartree/Bohr.
Within the :program:`Extf` module the keyword :kword:`LINEAR` is used. Note that :program:`Extf` needs to be called after the execution of :program:`ALASKA`, inside the loop. The options are:
(1) the atom number corresponding to the C atom
(2) the atom number corresponding to the N atom
(1.5) the force intensity
(0) to indicate a repulsive force.

.. extractfile:: ug/extf.input

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

  >> FOREACH ITER in (1 .. 1000)

  &SEWARD

  >> IF ( $ITER = 1 )

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

  &extf
   LINEAR
   1
   2
   1.5
   0

  &Dynamix
   VELVer
   DT= 10.0
   VELO= 3
   THER= 2
   TEMP=300

  >> END DO

.. xmldoc:: </MODULE>
