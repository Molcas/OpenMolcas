Quickstart Guide for |molcas|
=============================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

Introduction
------------

Running |molcas| requires a small number of operations.
This section of the manual, entitled "Quickstart Guide for |molcas|"
is aimed at those users who want to immediately
run a simple |molcas| calculation in order to become familiar with the program.
Basic hints are included which set the proper environment, build simple input files, run a calculation, and
subsequently extract information from the resulting output.

.. index::
   single: MOLCAS Environment

|molcas| Environment Setup
--------------------------

The environment variable (:variable:`MOLCAS`) and |molcas| driver (:file:`molcas`) must be defined in order to run |molcas|.
The :variable:`MOLCAS` environment variable points to the root directory of the |molcas| installation and
can be defined by the bash shell command ::

  export MOLCAS=/home/molcas/molcas.version

The location of the |molcas| driver is defined at installation time and is
typically located in :file:`/usr/local/bin` or :file:`$HOME/bin`.
Check to ensure that this directory is included in your path. Otherwise, the path can be extended
by the following command: ::

  export PATH=$PATH:$HOME/bin

In addition, the variable :variable:`MOLCAS_NPROCS` is needed to run |molcas| in parallel. This specifies the number of MPI processes that will be used.

It may be also convenient to define environment variables such as :variable:`WorkDir` which points to a directory for intermediate
files and :variable:`Project` to define the name of a project: ::

  export Project=MyMolecule

|molcas| will provide default values if they are not explicitly defined. For a discussion of other |molcas|
environment variables, please see :numref:`TUT:sec:environment`. All environment variables can
either be defined explicitly or entered in a shell script which can be subsequently executed.

.. index::
   single: Customization

Customization of |molcas| Execution
-----------------------------------

|molcas| has flexible control of organizing filenames and directories used during a calculation.
The default values used for customization can be altered either by shell variables or
a resource file :file:`molcasrc` which is more preferable. A command :command:`molcas setuprc`
provides guided help if to create such file.

The terminology used in this chapter:

* ``LOG``: the output and error files produced by |molcas|.

* ``ProjectName``: the Project name used for file naming.

* ``RUNFILE``: a file used in a calculation will be named as ``ProjectName``.Runfile,

* ``WorkDirName``: the WorkDir name used as the directory for temporary/binary files produced by |molcas|.

* ``Scratch``: the scratch disk area which provides a path to a parent directory for ``WorkDirName``\s.

  The :file:`WorkDir` variable used in the |molcas| manual is constructed as ``Scratch``/``WorkDirName``,

* ``CurrDir``: the submit directory where the |molcas| command was issued.

  Note, that in this tutorial, it is assumed that the input file is located in ``CurrDir``,

* ``OutputDir``: the output directory which is used for storage of extra output files, such as Orbital files and molden files.

It is quite important to understand, that if a user performs two consecutive runs of molcas, using the same
scratch area (:variable:`WorkDir`) and project name, |molcas| will try to reuse intermediate data, e.g.
integrals and orbitals, in order to make a restart of a calculation. This can save time, but can also be
can be dangerous if two consecutive calculations are not compatible.

Assuming that :file:`molcasrc` does not exist, and no environment is set, the command :command:`molcas inputfile`
will use the following defaults:

* ``LOG`` is printed to the screen,
* ``OutputDir`` and ``CurrDir`` are defined to be the same directory,
* ``ProjectName`` is s taken as the name of :file:`inputfile` by removing the suffix (before the last . (dot) character),
* ``Scratch`` is defined as :file:`/tmp/`,
* and ``WorkDirName`` is defined from the ``ProjectName`` plus a random suffix.

.. compound::

  For example, when a user issues the following commands: ::

    cd /home/joe/projects/water
    vi H2O.DFT.input
    molcas H2O.DFT.input

  the following files will be generated: ::

    /home/joe/projects/water/H2O.DFT.ScfOrb
    /home/joe/projects/water/H2O.DFT.scf.molden
    ...
    /tmp/H2O.DFT.15014/H2O.DFT.RunFile
    ...

If a flag :command:`-f` is used in a |molcas| command, ``LOG`` files will be stored in the ``CurrDir`` directory with a name
``ProjectName``.log and ``ProjectName``.err.

.. compound::

  ``ProjectName`` can either be set in a shell script running |molcas| or included directly into the |molcas| command: ::

    molcas Project=water H2O.DFT.input

  will change the default value for ``ProjectName`` to water.

If the :variable:`MOLCAS_WORKDIR` environment variable is set either as part of |molcas| command or is included in the :file:`molcasrc` file,
the name of WorkDir will NOT be random, but determined by the ``ProjectName``.

.. compound::

  Example: ::

    cd /home/joe/projects/water
    vi H2O.DFT.input
    molcas MOLCAS_WORKDIR=/tmp Project=water -f H2O.DFT.input

  will generate the following files: ::

    /home/joe/projects/water/water.log
    /home/joe/projects/water/water.ScfOrb
    ...
    /tmp/water/water.RunFile
    ...

For More options to control the behavior of |molcas|, run the command :command:`molcas setuprc` script.
The file :file:`molcasrc` can be used to set global preferences for the |molcas| package and/or to set user preferences.
The :file:`setuprc` script creates a :file:`molcasrc` file (:file:`HOME/.Molcas`) in a users home directory.

The following :file:`molcasrc` file for uses the :file:`/scratch` area as a parent for WorkDirs and
Project name generated for the the name of the input file,
then removes WorkDir before a calculation followed by subsequent
retains of this file when the calculation finished: ::

  # Version 1.0
  MOLCAS_MEM=256
  MOLCAS_WORKDIR=/scratch
  MOLCAS_NEW_WORKDIR=YES
  MOLCAS_KEEP_WORKDIR=YES
  MOLCAS_PROJECT=NAME

Once the :file:`molcasrc` is created, it is usually not necessary to use shell script or environment variables to run |molcas|.

|molcas| Command-Line Help System
---------------------------------

Just by typing :command:`molcas help` you get access to |molcas| Command-Line
Help System. There are different options:

* :command:`molcas help` produces a list of available programs and utilities.
* :command:`molcas help module` yields the list of keywords of the program :program:`module`.
* :command:`molcas help module keyword` offers the detailed description of the keyword.
* :command:`molcas help -t text` displays a list of keywords that contain the text word
  in their description.

.. index::
   single: EMIL commands
   single: MOLCAS input

Input Structure and EMIL Commands
---------------------------------

|molcas| has a modular program structure. The easiest way to run calculations
is to prepare an input file in which the different programs are executed
sequentially when the the module name (&module) is provided. The
keywords of module name then follow, with each entry on a separate line or
several entries on one line, separated by ;.
In addition to specific program module keywords, |molcas|
incorporates certain commands (See section on EMIL Commands.) that allow
operations such as looping over the modules, allowing partial execution,
changing variables, and substituting certain Unix commands.

Basic Examples
--------------

Simple Calculation on Water
...........................

.. compound::

  Start by preparing a file containing the cartesian coordinates of a water molecule. ::

    3
    Angstrom
     O       0.000000  0.000000  0.000000
     H       0.758602  0.000000  0.504284
     H       0.758602  0.000000 -0.504284

  which is given the name :file:`water.xyz`. In the same directory we prepare
  the input for the |molcas| run. We can name it :file:`water.input`.

In addition to using an editor to insert atomic coordinates into a file, a coordinate file can be obtained by using
a graphical interface program, for example, the :program:`LUSCUS` module as shown later in this guide. ::

  &GATEWAY
   coord=water.xyz
   basis=sto-3g
  &SEWARD
  &SCF

The :program:`GATEWAY` program module combines the molecular geometric of water
(In this case, from the external file, :file:`water.xyz`) and the basis set definition.
The :program:`SEWARD` program module then computes the integrals, and :program:`SCF` program modules
completer the calculation by computing the Hartree--Fock wave function.

.. For convenience just define: ::

    export Project=water

To run the calculation, the following command is used: ::

  molcas water.input -f

The file :file:`water.log` now contains output from the calculation, and the :file:`water.err`
includes any error messages. In the same directory, other files, including
:file:`water.scf.molden` or :file:`water.lus` (if the keyword :kword:`grid_it` is added at end of input file)
that help to analyze the results graphically with the external graphical viewer :program:`LUSCUS`
or :program:`Molden` program. Examples of their use are demonstrated below.

In the case of an open-shell calculation (UHF or UDFT), the :program:`SCF` program is again used.
Below, two examples are shown:

#. A UDFT calculation yielding an approximate doublet by setting the charge to +1, even if they are not pure spin functions: ::

     &GATEWAY
      coord=water.xyz
      basis=sto-3g
     &SEWARD
     &SCF
      charge=+1
      uhf; ksdft=b3lyp

#. A triplet state (using keyword :kword:`ZSPIn` to specify that there are two more :math:`\alpha` than :math:`\beta` electrons) states: ::

     &GATEWAY
      coord=water.xyz
      basis=sto-3g
     &SEWARD
     &SCF
      zspin=2
      uhf; ksdft=b3lyp

Geometry Optimization
.....................

In the next example, a DFT/B3LYP geometry optimization is performed on the
ground state of the water molecule. Notice that, after ``&gateway`` has defined
the coordinates and basis set definition, the EMIL commands :command:`>>> Do while`
and :command:`>>> EndDo` are employed to form a loop with the
:program:`seward`, :program:`SLAPAF`, and :program:`SCF` programs until convergence of geometry optimization is reached.
Program :program:`seward` computes the integrals in atomic basis, :program:`SCF` computes the DFT energy, and the program
:program:`SLAPAF` controls the geometry optimization and uses the module :program:`ALASKA` to compute the gradients
of the energy with respect to the degrees of freedom. :program:`SLAPAF` generates
the new geometry to continue the iterative structure optimization process and
checks to determine convergence parameters are satisfied notifying |molcas| and stopping the loop. ::

  &GATEWAY
   coord=water.xyz
   basis=ANO-S-MB
  >>> Do While
    &SEWARD
    &SCF
      ksdft=b3lyp
    &SLAPAF
  >>> EndDo

The above example illustrates the default situation of optimizing to a minimum geometry without
any further constraint. If other options are required such as determining a transition
state, obtaining a states crossing, or imposing a geometry constraint, specific input
should be added to program :program:`SLAPAF`.

.. figure:: acrolein.*
   :name: fig:ac1
   :align: center

   The acrolein molecule.

One of the most powerful aspects of |molcas| is the possibility of computing
excited states with multiconfigurational approaches. The next example demonstrates
a calculation of the five lowest singlet roots in a State-Average (SA) CASSCF calculation
using the :program:`RASSCF` program. It also illustrates the addition of the :program:`CASPT2` program
to determine dynamical correlation which provides accurate electronic energies at the CASPT2 level. The resulting
wave functions are used in the :program:`RASSI` module to calculate state-interaction properties such as oscillator strengths and other properties. ::

  &gateway
  Coord
   8
  Acrolein coordinates in Angstrom
   O     -1.808864   -0.137998    0.000000
   C      1.769114    0.136549    0.000000
   C      0.588145   -0.434423    0.000000
   C     -0.695203    0.361447    0.000000
   H     -0.548852    1.455362    0.000000
   H      0.477859   -1.512556    0.000000
   H      2.688665   -0.434186    0.000000
   H      1.880903    1.213924    0.000000
  Basis=ANO-S-MB
  Group=Nosym
  &SEWARD
  &RASSCF
    nactel  = 6 0 0
    inactive= 12
    ras2    = 5
    ciroot  = 5 5 1
  &CASPT2
    multistate=5 1 2 3 4 5
  &RASSI
    Nr of Job=1 5; 1 2 3 4 5
    EJob

Notice that the :kword:`Group` with the option :kword:`Nosym` has been used
to prevent :program:`GATEWAY` from identifying the symmetry of the molecule
(:math:`C_s` in this case). Otherwise, the input of the :program:`RASSCF` program
will have to change to incorporate the classification of the active space
into the corresponding symmetry species. Working with symmetry will be skipped at
this stage, although its use is very convenient in many cases.
A good strategy is to run only :program:`GATEWAY` and let the program guide you.

The :program:`RASSCF` input describes the active space employed, composed by
six active electrons distributed in five active orbitals. By indicating
twelve inactive orbitals (always doubly occupied), information
about the total number of electrons and the distribution of the orbitals is then complete.
Five roots will be obtained in the SA-CASSCF procedurei, and all them will
be computed at the CASPT2 level to obtain the transition energies at the higher
level of theory. Further, the :program:`RASSI` will compute the transition properties,
in particular, transition dipole moments and oscillator strengths.

GASSCF method
.............

In certain cases it is useful/necessary to enforce restrictions on electronic
excitations within the active space beyond the ones accessible by RASSCF.
These restrictions are meant to remove configurations that contribute only
marginally to the total wave function.
In |molcas| this is obtained by the GASSCF approach :cite:`gas2011`.
In GASSCF an arbitrary number of active spaces may be chosen.
All intra-space excitations are allowed (Full-CI in subspaces).
Constraints are imposed by user choice on inter-space excitations.
This method, like RASSCF, allows restrictions on the active space,
but they are more flexible than in RASSCF.
These restrictions are particularly useful when the cost of using the full CI
expansion of the active space is beyond reach.
These restrictions allow GASSCF to be applied to larger and more complex systems
at affordable cost.
Instead of a maximum number of holes in RAS1 and particles in RAS3, accumulated
minimum and maximum numbers of electrons are specified for GAS1, GAS1+GAS2,
GAS1+GAS2+GAS3, etc. in order to define the desired CI expansion.
The GAS scheme reduces to CAS or RAS when one or three spaces are chosen and
restrictions on electron excitations are adequately imposed.
When and how to use the GAS approach?
We consider three examples: (1) an organometallic material with separated metal
centers and orbitals not delocalized across the metal centers. One can include
the near degenerate orbitals of each center in its own GAS space.
This implies that one may choose as many GAS spaces as the number of
multiconfigurational centers. (2) Lanthanide or actinide metal compounds where
the :math:`f`-electrons require a MC treatment but they do not participate in bonding
neither mix with :math:`d` orbitals. In this case one can put the :math:`f` orbitals and their
electrons into one or more separated GAS spaces and not allow excitations
from and/or to other GAS spaces. (3) Molecules where each bond and its correlating
anti-bonding orbital could form a separate GAS space as in GVB approach.
Finally, if a wave function with a fixed number of holes in one or more
orbitals is desired, without interference of configurations where those
orbitals are fully occupied the GAS approach is the method of choice instead
of the RAS approach. There is no rigorous scheme to choose a GAS partitioning.
The right GAS strategy is system-specific. This makes the method versatile but
at the same time it is not a black box method.
An input example follow: ::

  &RASSCF
  nActEl
   6 0 0
  FROZen
  0 0 0 0 0 0 0 0
  INACTIVE
  2 0 0 0 2 0 0 0
  GASScf
  3
   1 0 0 0 1 0 0 0
  2 2
   0 1 0 0 0 1 0 0
  4 4
   0 0 1 0 0 0 1 0
  6 6
  DELEted
  0 0 0 0 0 0 0 0

In this example the entire active space counts six active electrons
and six active orbitals. These latter are partitioned in three GAS spaces
according to symmetry consideration and in the spirit of the GVB strategy.
Each subspace has a fixed number of electrons, *two*, and no interspace
excitations are allowed. This input shows clearly the difference
with the RAS approach.

Solvation Effects
.................

|molcas| incorporates the effects of the solvent using several models.
The most common is the cavity-based reaction-field Polarizable Continuum Model (PCM)
which is invoked by adding the keyword :kword:`RF-input` to the
:program:`SEWARD` code and is needed to compute the proper integrals. ::

  &GATEWAY
    coord=CH4.xyz
    Basis=ANO-S-MB
  &SEWARD
    RF-Input
     PCM-Model
     Solvent=Water
    End of RF-Input
  &RASSCF
    Nactel=8 0 0
    Inactive=1
    Ras2=8
  &CASPT2
    rfpert

The reaction field is computed in a self-consistent manner by the
:program:`SCF` or :program:`RASSCF` codes and added as a perturbation
to the Hamiltonian in the :program:`CASPT2` method with the keyword :kword:`RFPErt`.

Analyzing Results: Output Files and the :program:`LUSCUS` Program
-----------------------------------------------------------------

|molcas| provides a great deal of printed information in output files, and
the printing level is controlled by the environmental variable :variable:`MOLCAS_PRINT`.
By default this value is set to :kword:`two`, but can be modified by environmental variable :variable:`MOLCAS_PRINT`
Typical |molcas| output contains the program
header and input information, conditions of the calculation, the number of steps to achieve convergence, the energies and wave functions, and
final results, including in many cases the molecular orbital
coefficients as well as an analysis of the properties for the computed states.

.. For geometry optimizations, where many steps are required, there are different options to control
   how much output is generated. Three EMIL commands can be used:

   #. :command:`Set Output Screen` redirects the output to the screen;
   #. :command:`Set Output Over`, the default, skips the output of the intermediate steps and produces only output
      for the last iteration.
   #. :command:`Set Output File` places all output from
      each iteration in the :file:`$WorkDir` directory in files named
      :file:`Structure.iter.output`, where :command:`iter` is the number of
      the iteration.

      This is a convenient option to follow closely the convergence process. In this case, the :variable:`MOLCAS_PRINT` command must be set to :kword:`three`.

.. .. index::
      single: MING

   :program:MING:\: a Graphical Molcas Input Generator
   ...................................................

   |molcas| has incorporated a graphical self-guided tool to help the user to
   prepare the inputs and calculation flow named :program:`MING`. Provided that
   your system has all the graphical libraries and system utilities required
   for the code and this has been properly installed (try :command:`configure -ming`
   and read the installation guide if something fails), the :program:`MING`
   program is activated by the command :command:`molcas ming`.

   A window will then open in which the left panel contains three entrances.
   Clicking on :kword:`Templates` selected types of calculation are displayed
   in the right panel with the proper flow of |molcas| modules. Pressing on each
   of such boxes open new windows that enables to fill the input of the program.
   Most windows have a basic (default) and an advanced form. New modules or
   commands can be added to the flow by using the two other entrance in the
   left panel: :kword:`Commands`, including the EMIL commands that control
   the flow or add specific information, and :kword:`Modules`, including
   all |molcas| programs and utilities (see below).

   In the upper toolbar we find :command:`Settings`, to define the |molcas|
   environment, tools to :command:`Add`, :command:`Duplicate`, :command:`Delete`
   or :command:`Clear` new entrances, :command:`Preview` and :command:`Edit` the
   prepared input files, command :command:`Open` to retrieve previous input
   files, saving commands, and even commands to :command:`Submit` to send
   the calculation and command :command:`Watch` to inspect the output and
   error files.

   :program:`MING` can prepare most of calculations available in |molcas|.
   Even if you have a complex calculation the tool can be used to simplify
   and speed the basic aspects of the input.

.. index::
   single: LUSCUS

.. _TUT\:luscus:

LUSCUS: Grid and Geometry Visualization
.......................................

|molcas| developers have developed a graphical interface that can be used both
to create input for the |molcas| program and to analyze the results in
a graphical manner by visualizing molecular orbitals, density plots, and other output properties.

The first version of the code has the name GV (stands for Grid Viewer, or Geometry Visualization.
By an accident, the name also matches the nicknames of the main developers).
GV program uses a very limited set of graphic libraries, and thus has very primitive
user interface.

The next generation of GV program has the name LUSCUS. Luscus re-uses the code of GV,
and so GV users can use the same key combinations to operate with LUSCUS.
At the same time, LUSCUS provides a user-friendly interface, and contains many new
options, compared to GV.

LUSCUS can be obtained from http://luscus.sourceforge.net/, or
from http://www.molcas.org/LUSCUS.

LUSCUS can read the files only in one format: Luscus internal format (:file:`.lus`).
This format contains two sections: XYZ cartesian coordinates, and XML
formated data. It means that a standard XYZ file is a valid file in LUSCUS format.

Files with different formats, e.g. molden files, can be understood by LUSCUS
since they can be converted to LUSCUS format by a corresponding plug-in.
For instance, opening a file with the extension :file:`.molden`, LUSCUS automatically
runs a plug-in to convert a file from molden format to LUSCUS format. Saving a
LUSCUS file as a Molcas orbital file will automatically run a converter
from LUSCUS format to Orbital format.

* :command:`luscus xyz_file`: reads coordinates from a cartesian coordinate file.

  A molecule can be visualized and modified with the use of the
  left-button of the mouse and the keyboard. Below are some of the most
  useful commands.

  .. _tab\:luscus_geo:

  ======================== ===========================================================
  Left mouse click         Select atoms (if two, a bond is selected, if three
                           bond angle, if four a dihedral angle
  Left mouse + Shift click Mark/unmark atoms to/from the group
  Middle mouse/Space       Remove selection, or marking
  Insert key               Insert atom
  PageUp, PageDown         Alter type of selected atom or bond
  Delete/Supress key       Delete a selected atom
  +/-                      Change a value of selected bond/angle in steps
  Backspace                Undo last action
  Home                     Set selected atom to center of coordinates
  F8 key                   Find or apply symmetry
  ======================== ===========================================================

* :command:`luscus molden_file`: reads (check the comment about plug-in) from MOLDEN files such as
  :file:`wavefunction.molden`, :file:`freq.molden`, and :file:`geo.molden`.

  Note that |molcas| produces molden files with several extensions, so it is
  recommended to visualize these files by using :program:`Luscus`.

* :command:`luscus grid_file`: reads coordinates and densities and molecular orbitals from
  a binary :file:`grid_file`.

  This file is generated by :program:`GRID_IT` and, by default, placed in the :file:`$WorkDir` directory with the name
  :file:`$Project.lus`. The program allows displaying total densities, molecular orbitals, and charge density differences.

  If |molcas| and Luscus are installed locally, :program:`Luscus` can also be called from user input as shown in the following example: ::

    &GATEWAY
       coord = acrolein.xyz
       basis = ANO-L-MB
    &SEWARD
    &SCF
    &GRID_IT
    ALL

    * running external GUI program luscus

    ! luscus $Project.lus

    * User has to select active space and save GvOrb file!

    &RASSCF
    Fileorb=$CurrDir/$Project.GvOrb

  Note, that in the example above, the :program:`GRID_IT` program will generate a
  :file:`$Project.lus` file which :program:`LUSCUS` then uses, eliminating the need for defining
  :file:`$Project.lus` and allowing this file to be overwritten. :program:`rasscf` will
  read starting orbitals from the :file:`$Project.GvOrb` file.
