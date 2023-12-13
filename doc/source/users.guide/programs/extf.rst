o. index::
   single: Program; Extf
   single: Extf

.. _UG\:sec\:extf:

:program:`extf`
===============

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="EXTF">
            %%Description:
            <HELP>
            This module calculates the contribution of an external force that is acting on the system.
            </HELP>

This module calculates the contribution of an external force that is acting on the system. It applies the modification directly on the gradient and it needs to be called after the execution of :program:`ALASKA`, in an optimization or molecular dynamics calculation. The keyword :kword:`LINEAR` applies a constant linear force between two atoms :cite:`valentini2017optomechanical`.

.. _UG\:sec\:extf_inp:

Input
-----

General keywords
................

.. class:: keywordlist

:kword:`MODULE`
  Module of the force to apply, in nanonewton. If it's negative, the force is applied in opposite direction. See the other keywords for what is the direction of positive and negative forces. Note that this is the module of the total force, so for example, in the case of a force pair between two atoms, the force applied on each atom will be a factor of :math:`\sqrt{2}` smaller than this value.

  .. xmldoc:: <KEYWORD MODULE="EXTF" NAME="MODULE" APPEAR="Force module" KIND="REAL" LEVEL="BASIC">
              %%Keyword: MODUle <basic>
              <HELP>
              Module of the force to apply, in nanonewton. If it's negative, the force is applied in opposite direction. See the other keywords for what is the direction of positive and negative forces.
              </HELP>
              </KEYWORD>

:kword:`LINEAR`
  This keyword is followed by two integer values, specifying the atom numbers (following the numbering of the geometry) between which a force is applied along the vector joining them. A positive force (see the :kword:`MODULE` keyword) means an attractive (compression) force, a negative force is a repulsive (extension) force.

  .. xmldoc:: <KEYWORD MODULE="EXTF" NAME="LINEAR" APPEAR="Linear external force" KIND="INTS" SIZE="2" LEVEL="BASIC" EXCLUSIVE="BENDING,TORSIONAL">
              %%Keyword: LINEar <basic>
              <HELP>
              This keyword is followed by two integer values, specifying the atom numbers (following the numbering of the geometry) between which a force is applied along the vector joining them. A positive force (see the MODULE keyword) means an attractive (compression) force, a negative force is a repulsive (extension) force.
              </HELP>
              </KEYWORD>

:kword:`BENDING`
  This keyword is followed by three integer values, specifying the atom numbers (following the numbering of the geometry) between which a force is applied to open or close their planar angle. A positive force (see the :kword:`MODULE` keyword) tends to close the angle, a negative force opens it.

  .. xmldoc:: <KEYWORD MODULE="EXTF" NAME="BENDING" APPEAR="Bending external force" KIND="INTS" SIZE="3" LEVEL="BASIC" EXCLUSIVE="LINEAR,TORSIONAL">
              %%Keyword: BENDing <basic>
              <HELP>
              This keyword is followed by three integer values, specifying the atom numbers (following the numbering of the geometry) between which a force is applied to open or close their planar angle. A positive force (see the MODULE keyword) tends to close the angle, a negative force opens it.
              </HELP>
              </KEYWORD>

:kword:`TORSIONAL`
  This keyword is followed by four integer values, specifying the atom numbers (following the numbering of the geometry) between which a force is applied to open or close their dihedral angle. A positive force (see the :kword:`MODULE` keyword) tends to close positive dihedrals (i.e. towards less positive values), a negative force opens positive dihedrals (towards more positive values).

  .. xmldoc:: <KEYWORD MODULE="EXTF" NAME="TORSIONAL" APPEAR="Torsional external force" KIND="INTS" SIZE="3" LEVEL="BASIC" EXCLUSIVE="LINEAR,BENDING">
              %%Keyword: TORSional <basic>
              <HELP>
              This keyword is followed by three integer values, specifying the atom numbers (following the numbering of the geometry) between which a force is applied to open or close their planar angle. A positive force (see the MODULE keyword) tends to close the angle, a negative force opens it.
              </HELP>
              </KEYWORD>

:kword:`GAUSSIAN`
  This keyword modulates the applied force with a Gaussian time profile. It is followed by two real values, indicating the time at which the force is maximum (i.e. the value specified by :kword:`MODULE`) and a sigma value for the Gaussian decay.

  .. xmldoc:: <KEYWORD MODULE="EXTF" NAME="GAUSSIAN" APPEAR="Gaussian temporal profile" KIND="REALS" SIZE="2" LEVEL="BASIC">
              %%Keyword: GAUSsian <basic>
              <HELP>
              This keyword modulates the applied force with a Gaussian time profile. It is followed by two real values, indicating the time at which the force is maximum (i.e. the value specified by MODULE) and a sigma value for the Gaussian decay.
              </HELP>
              </KEYWORD>

Input examples
..............

The following input example is a semiclassical molecular dynamics with tully surface hop, where a linear force of about 2.9 nN is applied between atom 1 and atom 2. ::

  &Gateway
  coord=$Project.xyz
  basis=6-31G*
  group=nosym

  >> EXPORT MOLCAS_MAXITER=400
  >> DOWHILE

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
   1 2
   MODULE
   -4.1

  &Dynamix
   velver
   dt = 41.3
   velo = 1
   thermo = 0
  >>> End Do

This example shows an excited state CASSCF MD simulation
of a methaniminium cation using the Tully Surface Hop algorithm. In the simulation, the carbon and the nitrogen are pulled apart with a constant force of 1.5 nN (nanonewton) on each atom.
Within the :program:`Extf` module the keyword :kword:`LINEAR` is used. Note :program:`Extf` needs to be called after the execution of :program:`ALASKA`, inside the loop.

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

  &extf
   LINEAR
   1 2
   MODULE
   -2.12132

  &Dynamix
   VELVer
   DT= 10.0
   VELO= 3
   THER= 2
   TEMP=300

  >> END DO

.. xmldoc:: </MODULE>
