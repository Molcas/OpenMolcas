.. index::
   single: Program; MULA
   single: MULA

.. _sec\:mula:

:program:`Mula` |extramark|
===========================

.. warning::

   This program is not available in |openmolcas|

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: %%Description:
            This program computes intensities of vibrational
            transitions between electronic states.

The :program:`MULA` calculates intensities of vibrational
transitions between electronic states.

.. index::
   pair: Dependencies; MULA

.. _sec\:mula_dependencies:

Dependencies
------------

The :program:`MULA` program may need one or more UNSYM files produced
by the :program:`MCLR` program, depending on input options.

.. index::
   pair: Files; MULA

.. _sec\:mula_files:

Files
-----

Input files
...........

.. class:: filelist

:file:`UNSYM`
  Output file from the :program:`MCLR` program.

Output files
............

.. class:: filelist

:file:`plot.intensity`
  Contains data for plotting an artificial spectrum.

.. index::
   pair: Input; MULA

.. _sec\:mula_input:

Input
-----

The input for :program:`MULA` begins after the program name: ::

  &MULA

There are no compulsory keyword.

.. index::
   pair: Keywords; MULA

Keywords
........

.. class:: keywordlist

:kword:`TITLe`
  Followed by a single line, the title of the calculation.

  .. xmldoc:: %%Keyword: TITLe <basic>
              A single title line follows.

:kword:`FORCe`
  A force field will be given as input (or read from file), defining two
  oscillators for which individual vibrational levels and transition
  data will be computed.

  .. xmldoc:: %%Keyword: FORCe <basic>
              A force field will be given as input (or read from file).

:kword:`ATOMs`
  Followed by one line for each individual atom in the molecule.
  On each line is the label of the atom, consisting of an element symbol
  followed by a number. After the label, separated by one or more blanks,
  one can optionally give a mass number; else, a standard mass taken from
  the file data/atomic.data.
  After these lines is one single line with the keyword "END of atoms".

  .. xmldoc:: %%Keyword: ATOMs <basic>
              Followed by one line with an atom label for each individual atom
              in the molecule. A label consists of element name followed by a
              numeric label, optionally followed by a nuclear mass.0

:kword:`INTErnal`
  Specification of which internal coordinates that are to be used in the
  calculation. Each subsequent line has the form "BOND *a* *b*"
  or "ANGLE *a* *b* *c*" or
  or "TORSION *a* *b* *c* *d*" or
  or "OUTOFPL *a* *b* *c* *d*", for bond distances,
  valence angles, torsions (e.g. dihedral angles), and out-of-plane angles.
  Here, *a*...\ *d* stand for atom labels.
  After these lines follows one line with the keyword "END of internal".

  .. xmldoc:: %%Keyword: INTErnal <basic>
              Followed by lines of the form e.g. 'BOND C11 Br3', i.e. coordinate type
              and atom labels, Other choices are 'ANGLE a b c', 'TORSION a b c d'
              and 'OUTOFPL a b c d', where a--d are atom labels.

:kword:`MODEs`
  Selection of modes to be used in the intensity calculation. This is
  followed by a list of numbers, enumerating the vibrational modes to use.
  The modes are numbered sequentially in order of vibrational frequency.
  After this list follows one line with the keyword "END of modes".

  .. xmldoc:: %%Keyword: MODEs <basic>
              Selection of modes to be used in the intensity calculation.

:kword:`MXLEvels`
  Followed by one line with
  the maximum number of excitations in each of the two states.

  .. xmldoc:: %%Keyword: MXLEvels <basic>
              Followed by one line with max excitation level in the two states.

:kword:`VARIational`
  If this keyword is included, a variational calculation will be made,
  instead of using the default double harmonic approximation.

  .. xmldoc:: %%Keyword: VARIational <basic>
              Make a variational calculation, nor harmonic approximation.

:kword:`TRANsitions`
  Indicates the excitations to be printed in the output.
  Followed by the word FIRST on one line, then a list of numbers which
  are the number of phonons --- the excitation level --- to be distributed
  among the modes, defining the vibrational states of the first
  potential function (force field). Then similarly, after a line with
  the word SECOND, a list of excitation levels for the second state.

  .. xmldoc:: %%Keyword: TRANsitions <basic>
              Followed by the word FIRST, then a line with a list of
              the number of phonons to be distributed among the modes,
              for the first state, then similarly for second state.

:kword:`ENERgies`
  The electronic :math:`T_0` energies of the two states, each value is followed by
  either "eV" or "au".

  .. xmldoc:: %%Keyword: ENERgies <basic>
              The electronic T_0 energies of the two states, each value
              followed by "eV" or "au".

:kword:`GEOMetry`
  Geometry input. Followed by keywords FILE, CARTESIAN, or INTERNAL.
  If FILE, the geometry input is taken from UNSYM1 and UNSYM2.
  If CARTESIAN or INTERNAL, two sections follow, one headed by a line
  with the word FIRST, the other with the word SECOND. For the CARTESIAN
  case, the following lines list the atoms and coordinates. On each line
  is an atom label, and the three coordinates (:math:`x,y,z`). For the INTERNAL
  case, each line defines an internal coordinate in the same way as for
  keyword INTERNAL, and the value.

  .. xmldoc:: %%Keyword: GEOMetry <basic>
              Geometry input follows. Next line is FILE, CARTESIAN, or INTERNAL.
              Followed by FIRST, then coordinates, then SECOND, then coordinates.
              Format: See User's Guide.

:kword:`MXORder`
  Maximum order of transition dipole expansion. Next line is 0, if the
  transition dipole is constant, 1 if it is a linear function, etc.

  .. xmldoc:: %%Keyword: MXORder <basic>
              Next line is 0 for constant transition dipol, 1 for linear function.

:kword:`OSCStr`
  If this keyword is included, the oscillator strength, instead of the
  intensity, of the transitions will calculated.

  .. xmldoc:: %%Keyword: OSCStr <basic>
              Print oscillator strengths rather than intensities.

:kword:`BROAdplot`
  Gives the peaks in the spectrum plot an artificial halfwidth. The default
  lifetime is :math:`130\cdot 10^{-15}` s but this can be changed with keyword
  LIFEtime followd by the value.

  .. xmldoc:: %%Keyword: BROAdplot <basic>
              Enter life time (sec) to be used for lifetime broadening of
              artificial spectrum.

:kword:`NANOmeters`
  If this keyword is included, the plot file will be in nanometers.
  Default is in eV.

  .. xmldoc:: %%Keyword: NANOmeters <basic>
              If this keyword is included, the plot file will be in nanometers.
              Default is in eV.

:kword:`CM-1`
  If this keyword is included, the plot file will be in
  cm\ :math:`^{-1}`. Default is in eV.

  .. xmldoc:: %%Keyword: CM-1 <basic>
              If this keyword is included, the plot file will be in cm^-1.
              Default is in eV.

:kword:`PLOT`
  Enter the limits (in eV, cm\ :math:`^{-1}`, or in nm) for the plot file.

  .. xmldoc:: %%Keyword: PLOT <basic>
              Enter the limits (in eV, cm^-1, or in nm) for the plot file.

:kword:`VIBWrite`
  If this keyword is included, the vibrational levels of the two states will
  be printed in the output.

  .. xmldoc:: %%Keyword: VIBWrite <basic>
              Print vibrational levels in the output.

:kword:`VIBPlot`
  Two files, plot.modes1 and plot.modes2, will be generated, with pictures of
  the normal vibrational modes of the two electronic states.

  .. xmldoc:: %%Keyword: VIBPlot <basic>
              Generate files plot.modes1 and plot.modes2 picturing normal modes.

:kword:`HUGElog`
  This keyword will give a much more detailed output file.

  .. xmldoc:: %%Keyword: HUGElog <basic>
              Much more detailed output.

  .. :kword:`EXPANSION`
       This keyword indicates that the calculation will be aborted after
       the calculation of the expansion point.

:kword:`SCALe`
  Scales the Hessians, by multiplying with the scale factors following this keyword.

  .. xmldoc:: %%Keyword: SCALe <basic>
              Enter scale factors that will multiply the Hessians.

:kword:`DIPOles`
  Transition dipole data. If MXORDER=0 (see above), there follows a single line
  with :math:`x,y,z` components of the transition dipole moment. If MXORDER=1 there
  are an additional line for each cartesian coordinate of each atom, with the
  derivative of the transition dipole moment w.r.t. that nuclear coordinate.

  .. xmldoc:: %%Keyword: DIPOles <basic>
              Transition dipole data follows. A single line with x,y,z components,
              if MAXORDER=0. Else additional lines with gradient values.

:kword:`NONLinear`
  Specifies non-linear variable substitutions to be used in the definition of
  potential surfaces.

  .. xmldoc:: %%Keyword: NONLinear <advanced>
              Specifies non-linear variable substitutions in definition of potential functions.

:kword:`POLYnomial`
  Gives the different terms to be included in the fit of the polynomial
  to the energy data.

  .. xmldoc:: %%Keyword: POLYnomial <advanced>
              Specifies which polynomial terms that are used in modeling potential functions.

:kword:`DATA`
  Potential energy surface data.

  .. xmldoc:: %%Keyword: DATA <basic>
              Grid data follows. See manual for format.

Input example
.............

::

  &MULA

  Title
   Water molecule

  Atoms
   O1
   H2
   H3
  End Atoms

  Internal Coordinates
   Bond  O1 H2
   Bond  O1 H3
   Angle H3 O1 H2
  End Internal Coordinates

  MxLevels
    0  3

  Energies
   First
    0.0 eV
   Second
    3.78 eV

  Geometry
   Cartesian
   First
    O1     0.0000000000      0.0000000000     -0.5000000000
    H2     1.6000000000      0.0000000000      1.1000000000
    H3    -1.6000000000      0.0000000000      1.1000000000
   End
   Second
    O1     0.0000000000      0.0000000000     -0.4500000000
    H2     1.7000000000      0.0000000000      1.0000000000
    H3    -1.7000000000      0.0000000000      1.0000000000
   End

  ForceField
   First state
   Internal
    0.55 0.07 0.01
    0.07 0.55 0.01
    0.01 0.01 0.35
   Second state
   Internal
    0.50 0.03 0.01
    0.03 0.50 0.01
    0.01 0.01 0.25

  DIPOles
    0.20 0.20 1.20

  BroadPlot
  LifeTime
   10.0E-15

  NANO
  PlotWindow
   260 305

  End of input

::

  &MULA

  TITLe
   Benzene

  ATOMs
    C1
    C2
    C3
    C4
    C5
    C6
    H1
    H2
    H3
    H4
    H5
    H6
  End of Atoms

  GEOMetry
   file

  INTERNAL COORDINATES
   Bond    C1 C3
   Bond    C3 C5
   Bond    C5 C2
   Bond    C2 C6
   Bond    C6 C4
   Bond    C1 H1
   Bond    C2 H2
   Bond    C3 H3
   Bond    C4 H4
   Bond    C5 H5
   Bond    C6 H6
   Angle   C1 C3 C5
   Angle   C3 C5 C2
   Angle   C5 C2 C6
   Angle   C2 C6 C4
   Angle   H1 C1 C4
   Angle   H2 C2 C5
   Angle   H3 C3 C1
   Angle   H4 C4 C6
   Angle   H5 C5 C3
   Angle   H6 C6 C2
   Torsion C1 C3 C5 C2
   Torsion C3 C5 C2 C6
   Torsion C5 C2 C6 C4
   Torsion H1 C1 C4 C6
   Torsion H2 C2 C5 C3
   Torsion H3 C3 C1 C4
   Torsion H4 C4 C6 C2
   Torsion H5 C5 C3 C1
   Torsion H6 C6 C2 C5
  END INTERNAL COORDINATES

  VIBPLOT
   cyclic 4 1

  ENERGIES
   First
    0.0 eV
   Second
    4.51 eV

  MODES
   14 30 5 6 26 27 22 23 16 17 1 2 9 10
  END

  MXLE - MAXIMUM LEVEL of excitation (ground state - excited state)
    2 2

  MXOR - MAXIMUM ORDER in transition dipole.
    1

  OscStr

  Transitions
   First
    0
   Second
    0 1 2

  FORCEFIELD
   First
     file
   Second
     file

  DIPOLES
   file
