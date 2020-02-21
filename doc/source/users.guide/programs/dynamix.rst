.. index::
   single: Program; Dynamix
   single: Dynamix

.. _UG\:sec\:dynamix:

:program:`dynamix`
==================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="DYNAMIX">
            %%Description:
            <HELP>
            The DYNAMIX program allows to do molecular dynamics simulations using
            the velocity Verlet algorithm. It has also the capability to detect
            non-adiabatic transition using a surface hopping algorithm.
            </HELP>

The :program:`DYNAMIX` program performs molecular dynamics (MD)
simulations in |molcas|. Here the nuclei are moved according to the
classical Newton's equations which are solved numerically using the
velocity Verlet algorithm :cite:`swope:637`. The algorithm requires
coordinates, velocities and forces as input. :program:`DYNAMIX` can be
used with any electronic structure method in |molcas|. Also environmental
effects can be taken into account in the MD simulation: the solvent can be
considered implicitly using the reaction field keyword in :program:`GATEWAY`
or explicitly in hybrid QM/MM calculation which requires the :program:`ESPF`
program.

When multiple electronic states are involved in a MD simulation, a trajectory
surface hopping (TSH) algorithm allows non-adiabatic transitions between
different states. This TSH algorithm evaluates the change of the wavefunction
along the trajectory and induces a hop if certain criteria a met (for further
details read the :program:`RASSI` section). In the current implementation the
surface hopping algorithm can be used only with state averaged CASSCF
wavefunction. However, an extension for CASPT2 and other methods are in preparation.

The Tully algorithm is available in a separate module :program:`Surfacehop`.

.. _UG\:sec\:dynamix_dependencies:

Dependencies
------------

The coordinates and the forces are required by the :program:`DYNAMIX` program.
:program:`DYNAMIX` reads the initial coordinates from the :file:`RUNFILE` and
updates them in each iteration. In addition :program:`DYNAMIX` depends on the
:program:`ALASKA` program, since it generates forces.

.. _UG\:sec\:dynamix_files:

Files
-----

.. _UG\:sec\:dynamix_inp_files:

Input files
...........

.. class:: filelist

:file:`velocity.xyz`
  Contains the initial velocities of the MD simulation.

.. _UG\:sec\:dynamix_output_files:

Output files
............

.. class:: filelist

:file:`RUNFILE`
  Trajectory information such as current time, velocities, etc. are stored in this file.

:file:`md.xyz`
  The coordinates for each step of the MD trajectory are saved here.

:file:`md.energies`
  The potential, kinetic and total energies are written to this file. In case of multiple
  electronic states, the energies of all roots are saved.

.. _UG\:sec\:dynamix_inp:

Input
-----

This section describes the input syntax of :program:`DYNAMIX` in the |molcas| program
package. In general a MD simulation requires a :kword:`DoWhile` or :kword:`ForEach` loop which contains
several programs to compute the energy and :program:`ALASKA` for subsequent gradient
computation. The :program:`DYNAMIX` input begins with the program name,
and is followed by the only compulsory keyword :kword:`VELV` which specifies the
velocity Verlet algorithm: ::

  &DYNAMIX
  VELV

General keywords
................

.. class:: keywordlist

:kword:`VELVerlet`
  This keyword specifies the velocity Verlet algorithm :cite:`swope:637` to solve Newton's
  equations of motion. It's the only compulsory keyword in the program.

  .. xmldoc:: <KEYWORD MODULE="DYNAMIX" NAME="VELVER" APPEAR="Velocity Verlet algorithm" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: VELVerlet <basic>
              <HELP>
              Specifies the velocity Verlet algorithm for MD simulation.
              </HELP>
              </KEYWORD>

:kword:`DTime`
  Defines the :math:`\delta t` which is the time step in the MD simulation and which is
  used for the integration of Newton's equations of motion.
  The program expects the time to be given in floating point
  format and in atomic unit of time (1 a.u. of time = :math:`2.42\cdot10^{-17}` s). (Default = 10).

  .. xmldoc:: <KEYWORD MODULE="DYNAMIX" NAME="DT" APPEAR="Time step" KIND="REAL" LEVEL="BASIC" DEFAULT_VALUE="10.0" MIN_VALUE="0.0">
              %%Keyword: DTime <advanced>
              <HELP>
              Defines the time step of the MD simulation.
              </HELP>
              </KEYWORD>

:kword:`VELOcities`
  Specifies how the initial velocities are generated.
  This keyword is followed by an integer on the next line. The internal
  unit of the velocities is [bohr\ :math:`\cdot`\(a.u. of time)\ :math:`^{-1}`].

  .. container:: list

    **0** --- Zero velocities. (Default)

    **1** --- The velocities are read from the file :file:`$Project.velocity.xyz`
    in :file:`$WorkDir`. This file contains velocities in the xyz format given in the same
    order as the atoms in coordinate file. The unit of the velocities is [bohr\ :math:`\cdot`\(a.u. of time)\ :math:`^{-1}`].

    **2** --- This option allows to read in mass-weighted velocities from the
    file :file:`$Project.velocity.xyz` in [bohr\ :math:`\cdot\sqrt{\text{a.m.u.}}\cdot`\(a.u. of time)\ :math:`^{-1}`].

    **3** --- This option takes random velocities from a Maxwell--Boltzmann distribution, at
    a given temperature, assuming that every component of the velocity can be considered as an independent gaussian random variable.

  .. xmldoc:: <KEYWORD MODULE="DYNAMIX" NAME="VELO" APPEAR="Initial velocities" KIND="CHOICE" LIST="0: Zero,1: Read Cartesian,2: Read mass-weighted,3: Maxwell-Boltzmann" LEVEL="ADVANCED" DEFAULT_VALUE="0">
              %%Keyword: VELOcities <advanced>
              <HELP>
              Specifies the initial velocities.
              </HELP>
              </KEYWORD>

:kword:`THERmostat`
  Regulates the control of the temperature by scaling the velocities. The option
  is an integer given on the next line.

  .. container:: list

    **0** --- No velocity scaling. (Default)

    **1** --- The velocities are scaled in order to keep the total energy constant.

    **2** --- The velocities are scaled according to the Nosé--Hoover chain of thermostats algorithm, used to perform molecular symulation at
    constant temperature, resulting in statistics belonging to the canonical ensemble (NVT).

  .. xmldoc:: <KEYWORD MODULE="DYNAMIX" NAME="THER" APPEAR="Thermostat" KIND="CHOICE" LIST="0: No scaling,1: Constant energy,2: Nosé-Hoover" LEVEL="ADVANCED" DEFAULT_VALUE="0">
              %%Keyword: THERmostat <advanced>
              <HELP>
              Keyword for temperature control.
              </HELP>
              </KEYWORD>

:kword:`TEMPerature`
  Defines the numerical value of the temperature, which is used together with the Nosé--Hoover
  chain of thermostats to perform molecular dynamics at constant temperature. (Default = 298.15 K)

  .. xmldoc:: <KEYWORD MODULE="DYNAMIX" NAME="TEMP" APPEAR="Temperature of the simulation" KIND="REAL" LEVEL="ADVANCED" DEFAULT_VALUE="298.15" MIN_VALUE="0.0">
              %%Keyword: TEMPerature <advanced>
              <HELP>
              Keyword to state the temperature of the simulation.
              </HELP>
              </KEYWORD>

:kword:`HOP`
  Enables the trajectory surface hopping algorithm if the integer given in
  the next line is bigger than 0. The integer also specifies how many
  non-adiabatic transitions are allowed between electronic states.

  .. xmldoc:: <KEYWORD MODULE="DYNAMIX" NAME="HOP" APPEAR="Maximum number of surface hops" KIND="INT" LEVEL="BASIC" DEFAULT_VALUE="0" MIN_VALUE="0">
              %%Keyword: HOP <basic>
              <HELP>
              Specifies the maximum number of transitions between electronic states.
              </HELP>
              </KEYWORD>

:kword:`OUT`
  Enables dynamics in reduced dimensionality.
  This keyword is followed by an integer on the next line, which defines the number of nuclear coordinates to project out from the trajectory (default 0).
  The coordinates to project out are then read from the files `out.00X.xyz`, in the xyz format given in the same order as the atoms in coordinate file.

  .. xmldoc:: <KEYWORD MODULE="DYNAMIX" NAME="OUT" APPEAR="Number of coordinates to project out" KIND="INT" LEVEL="ADVANCED" DEFAULT_VALUE="0">
              %%Keyword: OUT <advanced>
              <HELP>
              Enables reduced dimensionality.
              </HELP>
              </KEYWORD>

:kword:`RESTART`
  This keyword allows to restart the trajectory at a given time.
  The time is given on the next line in atomic units.

  .. xmldoc:: <KEYWORD MODULE="DYNAMIX" NAME="RESTART" APPEAR="Restart the trajectory" KIND="REAL" LEVEL="ADVANCED">
              %%Keyword: RESTart <advanced>
              <HELP>
              Restarts the trajectory at a given time, which is specified on the next line.
              </HELP>
              </KEYWORD>

:kword:`H5RESTART`
  This keyword allows to restart a trajectory calculation from an HDF5 file.
  The name of the restart file is given on the next line.

  .. xmldoc:: <KEYWORD MODULE="DYNAMIX" NAME="H5RESTART" APPEAR="Restart the trajectory from a H5 file" KIND="STRING" LEVEL="ADVANCED">
              %%Keyword: H5REstart <advanced>
              <HELP>
              Restarts a trajectory calculation from an HDF5 file, whose name is given on the next line.
              </HELP>
              </KEYWORD>

Input examples
..............

The following example shows the input for an excited state CASSCF molecular dynamics
simulation of a methaniminium cation using the :program:`DYNAMIX` program. The :kword:`DoWhile` loop
allows 1000 steps with 10 a.u. of time step size which leads to a total duration of
242 fs. In the :program:`RASSCF` program the second root is selected for gradient
calculation using the keyword :kword:`MDRLXR`. This input assumes that the a
:file:`JOBIPH` file with orbitals is already given. In each iteration the :file:`JOBIPH`
is updated to achieve a fast convergence of the CASSCF wavefunction.
A Nosé--Hoover chain of thermostats, enabled with :kword:`THERmo`\=2, is used to
reproduce dynamics at constant temperature, where the initial velocities are
taken from a Maxwell--Boltzmann distribution at 300 K.

.. extractfile:: ug/DYNAMIX.input

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

  &ALASKA

  &DYNAMIX
   VELVer
   DT= 10.0
   VELO= 3
   THER= 2
   TEMP=300
   HOP= 1

  >> END DO

.. xmldoc:: <KEYWORD MODULE="DYNAMIX" NAME="VV_FIRST" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

.. xmldoc:: </MODULE>

Dynamixtools
------------

This tool can be found into the :file:`Tools/` folder and it will provide some general tools to manage molecular dynamics calculations. At the moment it can be used to generate intial conditions (geometries and momenta) following a Boltzmann distribution, based on a frequency calculation. It is working with a :file:`freq.molden` file (:file:`.h5` support coming soon...).

From the command prompt: ::

  $ python3 dynamixtools.py -h
  usage: dynamixtools.py [-h] [-s SEED] [-l LABEL] -i I [-b BOL] -t TEMP

  optional arguments:
    -h, --help            show this help message and exit
    -s SEED, --seed SEED  indicate the SEED to use for the generation of randoms
    -l LABEL, --label LABEL
                          label for your project (default is "geom")
    -i I, --input I       path of the frequency h5 or molden file
    -b BOL, --boltzmann BOL
                          number of initial condition following Boltzmann
                          distribution (default 1)
    -t TEMP, --temperature TEMP
                          temperature in kelvin for the initial conditions

Having a :file:`water.freq.molden` file, this is the command to generate 200 initial conditions using 3435432 as seed and a temperature of 300 kelvin: ::

  $ python3 dynamixtools.py -i water.freq.molden -t 300 -b 200 -s 3435432

