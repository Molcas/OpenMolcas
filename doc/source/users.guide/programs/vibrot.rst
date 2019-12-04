.. index::
   single: Program; VibRot
   single: VibRot

.. _UG\:sec\:vibrot:

:program:`vibrot`
=================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="VIBROT">
            %%Description:
            <HELP>
            This program computes the vibrational-rotational spectrum of a
            diatomic molecule. In addition, spectroscopic constants are computed.
            The program can also compute transition probabilities and lifetimes
            for excited states.
            </HELP>

The program :program:`VIBROT` is used to compute a vibration-rotation
spectrum for a diatomic molecule, using as input a potential
computed over a grid. The grid should be dense around equilibrium (recommended
spacing 0.05 au) and should extend to large distance (say 50 au) if
dissociation energies are computed.

The potential is fitted to an analytical form using cubic splines. The
ro-vibrational Schrödinger equation is then solved numerically
(using Numerov's method) for one vibrational state at a time and for a
number of rotational quantum numbers as specified by input. The
corresponding wave functions are stored on file
:file:`VIBWVS` for later use. The ro-vibrational energies
are analyzed in terms of spectroscopic constants. Weakly bound potentials can be
scaled for better numerical precision.

The program can also be fed with property functions, such as a dipole moment
curve. Matrix elements over the ro-vib wave functions for the property in
question are then computed. These results can be used to compute IR
intensities and vibrational averages of different properties.

:program:`VIBROT` can also be used to compute transition properties between
different electronic states. The program is then run twice to produce two files
of wave functions. These files are used as input in a third run, which will
then compute transition matrices for input properties. The main use is to
compute transition moments, oscillator strengths, and lifetimes for ro-vib
levels of electronically excited states. The asymptotic energy difference
between the two electronic states must be provided using the :kword:`ASYMptotic`
keyword.

.. index::
   pair: Dependencies; VibRot

.. _UG\:sec\:vibrot_dependencies:

Dependencies
------------

The :program:`VIBROT` is free-standing and does not depend on any
other program.

.. index::
   pair: Files; VibRot

.. _UG\:sec\:vibrot_files:

Files
-----

Input files
...........

The calculation of vibrational wave functions and spectroscopic
constants uses no input files (except for the standard input).
The calculation of transition properties uses
:file:`VIBWVS` files from two preceding
:program:`VIBROT` runs, redefined as
:file:`VIBWVS1` and
:file:`VIBWVS2`.

Output files
............

:program:`VIBROT` generates the file
:file:`VIBWVS` with vibrational wave functions for each :math:`v` and :math:`J` quantum
number, when run in the wave function mode. If requested :program:`VIBROT` can
also produce files :file:`VIBPLT` with the fitted potential and property
functions for later plotting.

.. index::
   pair: Input; VibRot

.. _UG\:sec\:vibrot_input:

Input
-----

This section describes the input to the :program:`VIBROT` program in the
|molcas| program system. The program name is ::

  &VIBROT

.. index::
   pair: Keywords; VibRot

Keywords
........

The first keyword to
:program:`VIBROT` is an indicator for the type of calculation
that is to be performed. Two possibilities exist:

.. class:: keywordlist

:kword:`ROVIbrational spectrum`
  :program:`VIBROT` will perform a vib-rot analysis and compute
  spectroscopic constants.

  .. xmldoc:: <SELECT MODULE="VIBROT" NAME="TYPE" APPEAR="Calculation type" CONTAINS="ROVIB,TRANS">

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="ROVIB" APPEAR="Start vib-rot analysis" KIND="SINGLE" LEVEL="BASIC" EXCLUSIVE="TRANS">
              %%Keyword: ROVIbrational <basic>
              <HELP>
              Perform a vib-rot analysis and compute spectroscopic constants.
              </HELP>
              </KEYWORD>

:kword:`TRANsition moments`
  :program:`VIBROT` will compute transition moment integrals
  using results from two previous calculations of the vib-rot wave
  functions. In this case the keyword :kword:`Observable` should be
  included, and it will be interpreted as the transition dipole moment.

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="TRANS" APPEAR="Compute transition moments" KIND="SINGLE" LEVEL="BASIC" EXCLUSIVE="ROVIB">
              %%Keyword: TRANsition <basic>
              <HELP>
              Compute transition moment integrals using previous vib-rot wave
              functions.
              </HELP>
              </KEYWORD>

  .. xmldoc:: </SELECT>

Note that only one of the above keywords can be used in a single
calculation. If none is given the program will only process the input
section.

After this first keyword follows a set of keywords, which are used to
specify the run. Most of them are optional.

The compulsory keywords are:

.. class:: keywordlist

:kword:`ATOMs`
  Gives the mass of the two atoms. Write mass number (an integer) and the
  chemical symbol Xx, in this order, for each of the two atoms in free format. If
  the mass numbers is zero for any atom, the mass of the most abundant isotope
  will be used. All isotope masses are stored in the program. You may introduce
  your own masses by giving a negative integer value to the mass number (one of
  them or both). The masses (in unified atomic mass units, or Da) are then read
  on the next (or next two) entry(ies). The isotopes of hydrogen can be given as
  H, D, or T.

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="ATOMS" APPEAR="The two atoms" KIND="CUSTOM" LEVEL="BASIC">
              %%Keyword: ATOMs <basic>
              <HELP>
              Read the mass number and chemical symbol of the atoms from the next line.
              If the mass number is zero the mass of the most abundant isotope will be
              used. Use a negative mass number to input the mass (in unified atomic mass
              units) in the next entry.
              </HELP>
              </KEYWORD>

:kword:`POTEntial`
  Gives the potential as an arbitrary number of lines. Each line
  contains a bond distance (in au) and an energy value (in au). A plot file of the
  potential is generated if the keyword
  :kword:`Plot` is added after the last energy input. One more entry should then follow
  with three numbers
  specifying the start and end value for the internuclear distance and
  the distance between adjacent plot points. This input must only be
  given together with the keyword :kword:`RoVibrational spectrum`.

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="POTE" APPEAR="Potential" KIND="CUSTOM" LEVEL="BASIC">
              <ALTERNATE KIND="STRING" />
              %%Keyword: POTEntial <basic>
              <HELP>
              Read the potential from a file (in au). Format: distance, value one pair on
              each line. Only together with vib-rot calculation.
              </HELP>
              </KEYWORD>

In addition you may want to specify some of the following optional
input:

.. class:: keywordlist

:kword:`TITLe`
  One single title line

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="TITLE" APPEAR="Title" KIND="STRING" LEVEL="BASIC">
              %%Keyword: TITLe <basic>
              <HELP>
              One single title line
              </HELP>
              </KEYWORD>

:kword:`GRID`
  The next entries give the number of grid points used in the numerical
  solution of the radial Schrödinger equation. The default value is
  199. The maximum value that can be used is 4999.

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="GRID" APPEAR="Numerical grid" KIND="INT" LEVEL="BASIC" DEFAULT_VALUE="199" MIN_VALUE="1" MAX_VALUE="4999">
              %%Keyword: GRID <basic>
              <HELP>
              Give the number of numerical grid points (default is 199, max is 4999).
              </HELP>
              </KEYWORD>

:kword:`RANGe`
  The next entry contains two distances Rmin and Rmax (in au) specifying
  the range in which the vibrational wave functions will be computed.
  The default values are 1.0 and 5.0 au. Note that these values most
  often have to be given as input since they vary considerably from one
  case to another. If the range specified is too small, the program will
  give a message informing the user that the vibrational wave function
  is large outside the integration range.

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="RANGE" APPEAR="Integration range" KIND="REALS" SIZE="2" LEVEL="BASIC" DEFAULT_VALUES="1.0,5.0">
              %%Keyword: RANGe <basic>
              <HELP>
              Give the range (Rmin-Rmax) in which the wave functions will be computed
              in atomic units. Default is 1.0-5.0 au.
              </HELP>
              </KEYWORD>

:kword:`VIBRational`
  The next entry specifies the number of vibrational quanta for which the
  wave functions and energies are computed. Default value is 3.

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="VIBR" APPEAR="Vibrational quanta" KIND="INT" LEVEL="BASIC" DEFAULT_VALUE="3" MIN_VALUE="1">
              %%Keyword: VIBRational <basic>
              <HELP>
              Specify the number of vibrational quanta (default is 3).
              </HELP>
              </KEYWORD>

:kword:`ROTAtional`
  The next entry specifies the range of rotational quantum numbers.
  Default values are 0 to 5. If the orbital angular momentum quantum
  number (:math:`m_\ell`) is non zero, the lower value will be adjusted to
  :math:`m_\ell` if the start value given in input is smaller than
  :math:`m_\ell`.

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="ROTA" APPEAR="Rotational quanta" KIND="INTS" SIZE="2" LEVEL="BASIC" DEFAULT_VALUES="0,5" MIN_VALUE="0">
              %%Keyword: ROTAtional <basic>
              <HELP>
              Specify the range of rotational quantum numbers (default is 0-5).
              </HELP>
              </KEYWORD>

:kword:`ORBItal`
  The next entry specifies the value of the orbital angular momentum
  (0, 1, 2, etc.). Default value is zero.

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="ORBI" APPEAR="Orbital angular momentum" KIND="INT" LEVEL="BASIC" DEFAULT_VALUE="0" MIN_VALUE="0">
              %%Keyword: ORBItal <basic>
              <HELP>
              Specify the orbital angular momentum:, 0, 1, 2, ... (default is 0).
              </HELP>
              </KEYWORD>

:kword:`SCALe`
  This keyword is used to scale the potential, such that the
  binding energy is 0.1 au. This leads to better precision in the numerical
  procedure and is strongly advised for weakly bound potentials.

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="SCALE" APPEAR="Scaled potential" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: SCALe <basic>
              <HELP>
              The potential will be scaled to a bond energy of 0.1 au.
              </HELP>
              </KEYWORD>

:kword:`NOSPectroscopic`
  Only the wave function analysis will be carried out but not the
  calculation of spectroscopic constants.

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="NOSP" APPEAR="No spectroscopic constants" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: NOSPectroscopic <advanced>
              <HELP>
              No calculation of spectroscopic constants.
              </HELP>
              </KEYWORD>

:kword:`OBSErvable`
  This keyword indicates the start of input for radial functions of observables
  other than the energy, for example the dipole moment function. The next line
  gives a title for this observable. An arbitrary number of input lines follows.
  Each line contains a distance and the corresponding value for the observable.
  As for the potential, this input can also end with the keyword :kword:`Plot`,
  to indicate that a file of the function for later plotting is to be constructed.
  The next line then contains the minimum and maximum R-values and the
  distance between adjacent points. When this input is given with the top keyword
  :kword:`RoVibrational spectrum` the program will compute matrix elements for
  vibrational wave functions of the current electronic state. Transition moment
  integrals are instead obtained when the top keyword is :kword:`Transition
  moments`. In the latter case the calculation becomes rather meaningless if
  this input is not provided. The program will then only compute the overlap
  integrals between the vibrational wave functions of the two states.
  The keyword :kword:`Observable` can be repeated up to ten times in a
  single run. All observables should be given in atomic units.

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="OBSE" APPEAR="Observable" KIND="CUSTOM" LEVEL="BASIC">
              %%Keyword: OBSErvable <basic>
              <HELP>
              Input for radial functions of observables (in au). The input is read from a
              file. The user is asked to read the users guide to learn how to construct
              this file.
              </HELP>
              </KEYWORD>

:kword:`TEMPerature`
  The next entry gives the temperature (in K) at which the vibrational
  averaging of observables will be computed. The default is 300 K.

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="TEMP" APPEAR="Temperature" KIND="REAL" LEVEL="ADVANCED" DEFAULT_VALUE="300.0" MIN_VALUE="0.0">
              %%Keyword: TEMPerature <advanced>
              <HELP>
              Temperature for vibrational averaging of observables (default is 300 K).
              </HELP>
              </KEYWORD>

:kword:`STEP`
  The next entry gives the starting value for the energy step used in
  the bracketing of the eigenvalues. The default value is 0.004 au
  (88 :math:`\text{cm}^{-1}`). This value must be smaller than the
  zero-point vibrational energy of the molecule.

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="STEP" APPEAR="Numerical step size" KIND="REAL" LEVEL="ADVANCED" DEFAULT_VALUE="0.004" MIN_VALUE="0.0">
              %%Keyword: STEP <advanced>
              <HELP>
              Give the starting value for the energy step used in bracketing eigenvalues.
              Should be smaller than the zero point energy (default is 0.004 au).
              </HELP>
              </KEYWORD>

:kword:`ASYMptotic`
  The next entry specifies the asymptotic energy difference between
  two potential curves in a calculation of transition matrix elements.
  The default value is zero atomic units.

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="ASYM" APPEAR="Asymptotic energy difference" KIND="REAL" LEVEL="BASIC" DEFAULT_VALUE="0.0">
              %%Keyword: ASYMptotic <basic>
              <HELP>
              Specify the asymptotic energy difference between two potential curves in a
              calculation of transition matrix elements (default is 0.00 au).
              </HELP>
              </KEYWORD>

:kword:`ALLRotational`
  By default, when the :kword:`Transition moments` keyword is given, only the
  transitions between the lowest rotational level in each vibrational state are
  computed. The keyword :kword:`AllRotational` specifies that the transitions
  between all the rotational levels are to be included. Note that this may result
  in a very large output file.

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="ALLR" APPEAR="All rotational levels" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: ALLRotational <advanced>
              <HELP>
              Include all rotational levels in a transition moments calculation.
              </HELP>
              </KEYWORD>

:kword:`PRWF`
  Requests the vibrational wave functions to be printed in the output file.

  .. xmldoc:: <KEYWORD MODULE="VIBROT" NAME="PRWF" APPEAR="Print wave functions" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: PRWF <advanced>
              <HELP>
              Requests the vibrational wave functions to be printed.
              </HELP>
              </KEYWORD>

Input example
.............

::

  &VIBROT
  RoVibrational spectrum
  Title = Vib-Rot spectrum for FeNi
  Atoms = 0 Fe 0 Ni
  Potential
   1.0 -0.516768
   1.1 -0.554562
   [...]
  Plot  = 1.0 10.0 0.1
  Grid  = 150
  Range = 1.0 10.0
  Vibrations = 10
  Rotations  = 2 10
  Orbital    = 2
  Observable
   Dipole Moment
   1.0 0.102354
   1.1 0.112898
   [...]
  Plot  = 1.0 10.0 0.1

**Comments**: The vibrational-rotation spectrum for :math:`\ce{FeNi}`
will be computed using the potential curve given in input. The 10
lowest vibrational levels will be obtained and for each level the
rotational states in the range :math:`J`\=2 to 10. The vib-rot matrix elements
of the dipole function will also be computed. A plot file of the
potential and the dipole function will be generated. The masses for
the most abundant isotopes of :math:`\ce{Fe}` and :math:`\ce{Ni}` will be selected.

.. xmldoc:: </MODULE>
