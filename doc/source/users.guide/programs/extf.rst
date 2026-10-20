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

The :program:`EXTF` module modifies the molecular gradient computed by :program:`ALASKA` by adding a contribution of an external mechanical force :cite:`valentini2017optomechanical`. :program:`EXTF` must be called after :program:`ALASKA`, in a geometry optimization or inside a molecular dynamics loop. In molecular dynamics simulations, :program:`EXTF` should be placed immediately before :program:`DYNAMIX`.

The general input structure is: ::

  &ALASKA
  &EXTF
    force_keyword
    force_parameters
  &DYNAMIX

Standard forces included are :kword:`LINEAR`, :kword:`BENDING`, and :kword:`TORSIONAL`, which apply external forces along a bond, an angle or a dihedral coordinate, respectively. The force magnitude is specified with the :kword:`MODULE` keyword and is given in nanonewton.
The keyword :kword:`WALL` is intended to mimic an external pressure acting on the molecular structure typically simulating ball-milling :cite:`DeArmas2025`. In this case, the :kword:`MODULE` keyword specifies the pressure in gigapascal.

.. _UG\:sec\:extf_inp:

Input
-----

General keywords
................

.. class:: keywordlist

:kword:`MODULE`
  Module of the force to apply, in nanonewton. If it's negative, the force is applied in opposite direction. See the other keywords for what is the direction of positive and negative forces. Note that this is the module of the total force, so for example, in the case of a force pair between two atoms, the force applied on each atom will be a factor of :math:`\sqrt{2}` smaller than this value.
  When using :kword:`WALL`, this is instead the pressure, in gigapascal.

  .. xmldoc:: <KEYWORD MODULE="EXTF" NAME="MODULE" APPEAR="Force module" KIND="REAL" LEVEL="BASIC">
              %%Keyword: MODUle <basic>
              <HELP>
              Module of the force to apply, in nanonewton. If it's negative, the force is applied in opposite direction. See the other keywords for what is the direction of positive and negative forces.
              With WALL, this is the pressure, in gigapascal.
              </HELP>
              </KEYWORD>

:kword:`LINEAR`
  This keyword is followed by two integer values, specifying the atom numbers (following the numbering of the geometry) between which a force is applied along the vector joining them. A positive force (see the :kword:`MODULE` keyword) means an attractive (compression) force, a negative force is a repulsive (extension) force.

  .. xmldoc:: <KEYWORD MODULE="EXTF" NAME="LINEAR" APPEAR="Linear external force" KIND="INTS" SIZE="2" LEVEL="BASIC" EXCLUSIVE="BENDING,TORSIONAL,WALL">
              %%Keyword: LINEar <basic>
              <HELP>
              This keyword is followed by two integer values, specifying the atom numbers (following the numbering of the geometry) between which a force is applied along the vector joining them. A positive force (see the MODULE keyword) means an attractive (compression) force, a negative force is a repulsive (extension) force.
              </HELP>
              </KEYWORD>

:kword:`BENDING`
  This keyword is followed by three integer values, specifying the atom numbers (following the numbering of the geometry) between which a force is applied to open or close their planar angle. A positive force (see the :kword:`MODULE` keyword) tends to close the angle, a negative force opens it.

  .. xmldoc:: <KEYWORD MODULE="EXTF" NAME="BENDING" APPEAR="Bending external force" KIND="INTS" SIZE="3" LEVEL="BASIC" EXCLUSIVE="LINEAR,TORSIONAL,WALL">
              %%Keyword: BENDing <basic>
              <HELP>
              This keyword is followed by three integer values, specifying the atom numbers (following the numbering of the geometry) between which a force is applied to open or close their planar angle. A positive force (see the MODULE keyword) tends to close the angle, a negative force opens it.
              </HELP>
              </KEYWORD>

:kword:`TORSIONAL`
  This keyword is followed by four integer values, specifying the atom numbers (following the numbering of the geometry) between which a force is applied to open or close their dihedral angle. A positive force (see the :kword:`MODULE` keyword) tends to close positive dihedrals (i.e. towards less positive values), a negative force opens positive dihedrals (towards more positive values).

  .. xmldoc:: <KEYWORD MODULE="EXTF" NAME="TORSIONAL" APPEAR="Torsional external force" KIND="INTS" SIZE="3" LEVEL="BASIC" EXCLUSIVE="LINEAR,BENDING,WALL">
              %%Keyword: TORSional <basic>
              <HELP>
              This keyword is followed by three integer values, specifying the atom numbers (following the numbering of the geometry) between which a force is applied to open or close their planar angle. A positive force (see the MODULE keyword) tends to close the angle, a negative force opens it.
              </HELP>
              </KEYWORD>

:kword:`WALL`
  This keyword applies a wall-type compression force to the system. If three real values are provided, these specify the normal vector. If only two real values are provided, these are the Euler :math:`\theta` and :math:`\phi` angles of the normal vector, in radians, or if the word ``DEGREE`` is added, in degrees. Note that if :kword:`WALL` is used, the :kword:`MODULE` keyword specifies the pressure in gigapascal. :cite:`DeArmas2025`

  .. xmldoc:: <KEYWORD MODULE="EXTF" NAME="WALL" APPEAR="Wall-type compression" KIND="REALS" SIZE="3" LEVEL="BASIC" EXCLUSIVE="LINEAR,BENDING,TORSIONAL">
              <ALTERNATE KIND="STRING" />
              %%Keyword: WALL <basic>
              <HELP>
              This keyword applies a wall-type compression force to the system.
              The orientation of the wall can be given in two different forms. If three real values are provided, these specify the normal vector. If only two real values are provided, these are the Euler theta and phi angles of the normal vector, in radians, or if the word DEGREE is added, in degrees.
              </HELP>
              </KEYWORD>

:kword:`VDWR`
  Optional keyword to be used together with :kword:`WALL`. It activates the use of van der Waals radii in the construction of the effective molecular cross section to determine the forces given the external pressure.

  .. xmldoc:: <KEYWORD MODULE="EXTF" NAME="VDWR" APPEAR="Van der Waals radii" KIND="SINGLE" LEVEL="BASIC" REQUIRE="WALL">
              %%Keyword: VDWR <basic>
              <HELP>
              Optional keyword to be used together with WALL. It activates the use of van der Waals radii in the construction of the effective molecular cross section to determine the forces given the external pressure.
              </HELP>
              </KEYWORD>

:kword:`ROTA`
  Optional keyword to be used together with :kword:`WALL`. It disables the removal of the net torque, the system will be allowed to rotate.

  .. xmldoc:: <KEYWORD MODULE="EXTF" NAME="ROTA" APPEAR="Rotation allowed" KIND="SINGLE" LEVEL="ADVANCED" REQUIRE="WALL">
              %%Keyword: ROTA <advanced>
              <HELP>
              Optional keyword to be used together with WALL. It disables the removal of the net torque, the system will be allowed to rotate.
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
   angstrom
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

A minimal molecular dynamics example using :kword:`WALL` is: ::

  &GATEWAY
    Coord = molecule.xyz
    Basis = STO-3G
    Group = NoSym

  > export MOLCAS_MAXITER = 300
  > do while

  &SEWARD
  &SCF
  &ALASKA
  &EXTF
    Wall = 90.0 0.0 degree
    Module = 0.1
    VDWR
  &DYNAMIX
    VelVer
    dt = 10.0
    Velocities = 3
    Thermostat = 0

  > end do

.. xmldoc:: </MODULE>
