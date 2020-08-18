.. index::
   single: Program; Geo
   single: Geo

.. _UG\:sec\:geo:

:program:`geo` |extramark|
==========================

.. warning::

   This program is not available in |openmolcas|

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="GEO" ONSCREEN="HF,HF-GEOMETRY,CAS,CAS-GEOMETRY">
            %%Description:
            <HELP>
            The GEO module handles geometry optimization in
            constrained internal coordinates.
            </HELP>

.. _UG\:sec\:geo_description:

Description
-----------

The :program:`geo` module handles geometry optimization in constrained
internal coordinates :cite:`Vysotskiy:13`. The module is called automatically when the
emil command :kword:`> DO GEO` is used and should never be called
explicitly by the user. All input relevant for these type of geometry
optimization should be supplied in :program:`gateway` and all the relevant
keywords are described in detail in the manual section for :program:`gateway`.

The purpose of the :program:`geo` module is to perform geometry optimization of
several molecular fragments by keeping the fragments rigid and only optimize
their relative position. The xyz-file of each fragment is supplied to the
program by a separate :kword:`coord` input in :program:`gateway`. Internal
coordinates for the whole complex is then constructed and stored into the
file :file:`$project.zmt` of the
directory :file:`$GeoDir` which is a communication directory set up and
used by the :program:`geo` module. The internal coordinates are chosen so that
only a maximum of six coordinates link each fragment. All coordinates within
fragments are frozen and the optimization is only carried out for these
linking coordinates. The geometry optimization is performed using a numerical
gradient (and hessian if needed) and could in principle be used together
with any energy
relaxation method that is implemented into molcas. (Currently it only works
with :program:`scf`, :program:`mp2`, :program:`rasscf`, :program:`caspt2`,
:program:`chcc`, :program:`cht3`, and :program:`ccsdt` but it is very simply to extend
it to other modules.)

The module is intended to use in cases where one knows
the geometry of each fragment and do not expect it to change much during the
optimization.

.. compound::

  Note that it is often advantageous to run a GEO job in parallel. The number of
  processors (cores) required can be easily evaluated by the following formula:cite:`Vysotskiy:13`:

  .. math:: N_{\text{procs}} = 2 + \frac{N_{\text{var}}\left(N_{\text{var}}+3\right)}{2},

  where :math:`N_{\text{var}}` is the total number of active coordinates. In the simplest
  case of the :math:`N_{\text{frag}}` rigid polyatomic fragments, when only relative positions of fragments are
  going to be optimized, the :math:`N_{\text{var}}` parameter reads:

  .. math:: N_{\text{var}} = 6 N_{\text{frag}}

  In particular, for constrained two-fragment geometry optimization one needs 29 processes
  to run GEO job in a fully parallel manner.

Creating the z-matrix
.....................

The z-matrix is created by choosing the first atom in the xyz-file and put that
in top of the z-matrix file (the example is a methane dimer): ::

  H

The next atom is chosen as the closest atom not already in the z-matrix and in
the same fragment
as the first. The distance to the closest neighbor within the z-matrix is
calculated and written into the zmatrix: ::

  C
  H   1.104408  0                         1

The number zero means "do not optimize" and is appended to all coordinates
within fragments and the 1 at the right hand side keeps track of which
atom that were the closest neighbour. The next molecule is again chosen
as the one not in the z-matrix but in the fragment. The distance between
the third atom and its closest neighbour and the angle between the third
atom, the closest atom and the second closest atom with the second closest
in the middle is added. ::

  C
  H     1.104408 0                                 1
  H     1.104438 0   109.278819 0                  1   2

For the fourth atom a dihedral angle is added and after that all remaining atoms
only have coordinates related to three neighbors. When all the molecule in the
fragment has entered the z-matrix the closest atom in another fragment is taken
as the next atom. Neighbors to atoms in that fragment is chosen from the same
fragment if possible and then from the previous fragment. When two atoms from
the second fragment has been added the z-matrix looks like this: ::

  C
  H     1.104408 0                                 1
  H     1.104438 0   109.278819 0                  1   2
  H     1.104456 0   109.279688 0  -119.533195 0   1   2   3
  H     1.104831 0   109.663377 0   120.232913 0   1   2   3
  H     3.550214 1    75.534230 1  -133.506752 1   3   4   1
  C     1.104404 0    78.859588 1  -146.087839 1   7   3   3

Note that coordinates linking fragments have a 1 after and will be optimized
and that the last carbon has one neighbor within its own fragment and two
within the previous fragment.

There are two special exceptions to the rules described above. Firstly, if
it is possible neighboring atoms are chosen from non-hydrogen atoms. (When
linking an atom to neighbors, not when choosing the next line in the z-matrix.)
Secondly, dihedral angles are rejected if the angles between the three first or
the three last atoms are smaller than 3 degrees since very small angles makes
the dihedral ill-defined.

.. _UG\:sec\:geo_dependencies:

Dependencies
------------

The :program:`GEO` program requires that :program:`gateway` have been run with
either the keyword :kword:`geo` or :kword:`hyper` to setup internal
communication files.

.. _UG\:sec\:geo_files:

Files
-----

Input files
...........

Apart from the standard file :program:`Geo` will use the following input files.

.. class:: filelist

:file:`RUNFILE`
  File for communication of auxiliary information.

:program:`Geo` will also use internal communication files in the directory
:file:`$GeoDir` described in more detail in the next section.

Geo communication files
.......................

When the emil command :kword:`> DO GEO` or the keywords zonl or zcon
are used a directory
:file:`$GeoDir` is created (by default in the input-directory:
:file:`$CurrDir/$Project.Geo`).
This directory
is used to store files related to geometry optimization and z-matrix generation.

.. class:: filelist

:file:`$project.zmt`
  The file with the z-matrix as described above.

:file:`general.info`
  A file used for storing general info about the geometry optimization, it is human readable with labels.
  The file is automatically setup by the program.

  example: ::

    Number of iterations:              3
    Number of atoms:                   8
    Number of internal coordinates:    6
    Internal Coordinates:
        2.159252   99.560213  123.714490   99.612396 -179.885031 -123.791319
    Displacement parameters:
        0.150000    2.500000    2.500000    2.500000    2.500000    2.500000
    Coordinate types:
    badadd

  Most of the lines are self-explanatory, coordinate-type is one character for each
  internal coordinate to optimize with b=bond, a=angle and d=dihedral, displacement
  parameters is the coordinates defined by hyper.

:file:`disp????.info`
  A file that contains all displaced coordinates. A new instance of the file is created automatically for each geometry step.

  example: ::

    disp0001.xyz 2.16485  99.52171 123.71483  99.57425 -179.88493 -123.79124
    disp0002.xyz 2.17548 101.22524 124.91941 100.50733 -179.12307 -123.14735
    disp0003.xyz 2.16485 101.22952 124.92244 100.50968 -179.12116 -123.14573
    disp0004.xyz 2.17145  98.46360 125.21123 100.73337 -178.93851 -122.99136

:file:`$project.disp????.energy`
  A one line file that simply state the current energy for the displacement.
  If several energies are calculated this will first contain a scf energy
  and is then updated with the mp2 value when that calculation is finished.
  The file is written by a small addition to the add-info files and currently
  collects energies from :program:`scf`, :program:`mp2`, :program:`rasscf`, :program:`caspt2`,
  :program:`chcc`, :program:`cht3`, and :program:`ccsdt`. If
  it should be used with new energy relaxation methods this must be added manually.

:file:`$project.final????.xyz`
  An xyz-file with coordinates for all fragments after XXXX geometry optimizations.
  This is just for the benefit of the user and should probably be replaced with the
  same output as created by slapaf, opt.xyz-file for optimized geometry and molden-file
  for history. Some internal history would still be needed for building more advanced
  geometry optimization algorithms and convergence criteria though.

:file:`$project.geo.molden`
  A molden file with information about the geometry optimization.
  The file could be browsed in molden or LUSCUS. The last
  energy-value is set to zero since the file cannot be created
  after the last energy calculation and need to be inserted by the user.

.. _UG\:sec\:geo_input:

Input
-----

The general input structure of a geo-calculation looks like this: ::

  &Gateway
  [keywords to modify fragment position (frgm,origin)]
  coord=fragment1.xyz
  ...
  coord=fragmentN.xyz
  [regular :program:`gateway` keywords + keywords to modify :program:`geo`]

  >> DO GEO
  &Seward
  [energy relaxation methods with any of their keywords]
  >> END DO

Both the keywords used to translate and rotate xyz-files (:kword:`frgm`
and :kword:`origin`) and the keyword to modify the behaviour of the optimization
(:kword:`hyper`, and :kword:`oldz`) is
described in more details in the :program:`gateway` section of the manual.

Here is an example input to calculate the relative orientation of two methane
molecules:

.. extractfile:: ug/GEO.input

  >> export GeoDir=$CurrDir/$Project.GEO
  &Gateway
  Coord=$MOLCAS/Coord/Methane1.xyz
  Coord=$MOLCAS/Coord/Methane2.xyz
  Group=c1
  basis=aug-cc-pVDZ
  hyper
  0.15 2.5 2.5

  >> Do Geo

  &Seward
  CHOLESKY HIGH

  &SCF

  &MP2

  >> End Do

In order to run the job above in a task-farm parallel mode,
one just needs to slightly modify the original input:

.. extractfile:: ug/GEO.TASKFARM.input

  *activating an TaskFarm's interface in GEO/HYPER
  >> export MOLCAS_TASK_INPUT=NEW
  *specifying number of cores/cpus available for parallel execution
  >> export MOLCAS_NPROCS=2
  >> export GeoDir=$CurrDir/$Project.GEO
  &Gateway
  Coord=$MOLCAS/Coord/Methane1.xyz
  Coord=$MOLCAS/Coord/Methane2.xyz
  Group=c1
  basis=aug-cc-pVDZ
  hyper
  0.15 2.5 2.5

  >> Do Geo

  &Seward
  CHOLESKY HIGH

  &SCF

  &MP2

  >> End Do

  *setting the TaskMode's execution model to MPP
  *each energy/displacement is computed in paralel
  *by using 2 cores
  >>> EXPORT MOLCAS_TASKMODE=1
  *actual execution of TaskFarm via EMIL
  >>> UNIX SERIAL $MOLCAS/sbin/taskfarm

.. xmldoc:: </MODULE>
