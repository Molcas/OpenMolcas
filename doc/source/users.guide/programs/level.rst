.. index::
   single: Program; LEVEL
   single: LEVEL

.. _UG\:sec\:level:

:program:`level`
=================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="LEVEL">
            %%Description:
            <HELP>
            This program computes the vibrational-rotational spectrum of a
            diatomic molecule. In addition, spectroscopic constants are computed.
            The program can also compute expectation values and Franck-Condon 
            factors.
            </HELP>

The program :program:`LEVEL` is used to compute a vibration-rotation
spectrum for a diatomic molecule, using as input a potential
that is computed over a grid, or an analytic potential with its parameters 
specified. The grid should be dense around equilibrium (recommended
spacing 0.05 au) and should extend to a large distance (say 50 au) if
dissociation energies are computed.

The ro-vibrational Schrödinger equation is solved numerically
(using Numerov's method).  The ro-vibrational energies
are analyzed in terms of spectroscopic constants. 

.. index::
   pair: Dependencies; LEVEL

.. _UG\:sec\:level_dependencies:

Dependencies
------------

The :program:`LEVEL` is free-standing and does not depend on any
other program.

.. index::
   pair: Files; LEVEL

.. _UG\:sec\:level_files:

Files
-----

Input files
...........

The calculation of vibrational wavefunctions and spectroscopic
constants uses no input files (except for the standard input).

Output files
............

:program:`LEVEL` generates a standard output file which ends 
with a summary of all levels found.

.. index::
   pair: Input; LEVEL

.. _UG\:sec\:level_input:

Input
-----

This section describes the input to the :program:`LEVEL` program in the
|molcas| program system. The program name is ::

  &LEVEL

.. index::
   pair: Keywords; LEVEL

Keywords
........

The compulsory keywords are:

.. class:: keywordlist

:kword:`IAN1`
  Integer Atomic Number of atom 1

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="IAN1" KIND="INT" LEVEL="BASIC">
              %%Keyword: IAN1 <basic>
              <HELP>
              Read the integer atomic number of atom 1.
              </HELP>
              </KEYWORD>

:kword:`IMN1`
  Integer Mass Number of atom 1

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="IMN1" KIND="INT" LEVEL="BASIC">
              %%Keyword: IMN1 <basic>
              <HELP>
              Read the integer mass number of atom 1.
              </HELP>
              </KEYWORD>

:kword:`IAN2`
  Integer Atomic Number of atom 2

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="IAN2" KIND="INT" LEVEL="BASIC">
              %%Keyword: IAN2 <basic>
              <HELP>
              Read the integer atomic number of atom 2.
              </HELP>
              </KEYWORD>

:kword:`IMN2`
  Integer Mass Number of atom 2

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="IMN2" KIND="INT" LEVEL="BASIC">
              %%Keyword: IMN2 <basic>
              <HELP>
              Read the integer mass number of atom 2.
              </HELP>
              </KEYWORD>

:kword:`CHARge`
  Charge of molecule

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="CHARGE" KIND="INT" LEVEL="BASIC">
              %%Keyword: CHARge <basic>
              <HELP>
              Read the integer charge of the molecule.
              </HELP>
              </KEYWORD>

:kword:`NUMPot`
  Number of potentials

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="NUMPOT" KIND="INT" LEVEL="UNDOCUMENTED">
              %%Keyword: NUMPot <undocumented>
              <HELP>
              Number of potentials (1 for a single potential, 2 for two potentials and 
              calculation of matrix elements coupling their levels.
              </HELP>
              </KEYWORD>

:kword:`RH`
  Step size, :math:`\Delta R` for the numerical solution of the differential equation. Calculations should be done with smaller and smaller values of this variable (with all other variables kept the same) until convergence with respect to this variable is achieved. 

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="RH" KIND="REAL" LEVEL="BASIC">
              %%Keyword: RH <basic>
              <HELP>
              Read the real number value for the step size used for the numerical
              solution of the differential equation.
              </HELP>
              </KEYWORD>

:kword:`RMIN`
  Minimum value of :math:`R` for the numerical solution of the differential equation. Calculations should be done with smaller and smaller values of this variable (with all other variables kept the same) until convergence with respect to this variable is achieved.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="RMIN" KIND="REAL" LEVEL="BASIC">
              %%Keyword: RMIN <basic>
              <HELP>
              Read the real number value for the minimum value of R for the 
              numerical solution of the differential equation.
              </HELP>
              </KEYWORD>

:kword:`pRV`
  The :math:`p` value (power) for the "radial variable" used for numerically solving the differential equation

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="pRV" KIND="REAL" LEVEL="BASIC">
              %%Keyword: pRV <basic>
              <HELP>
              Read the power p for the radial "variable" used for
              numerically solving the differential equation.
              </HELP>
              </KEYWORD>

:kword:`aRV`
  The real number :math:`R` value around which the "raidial variable" used for numerically solving the differential equation, is centered.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="aRV" KIND="REAL" LEVEL="BASIC">
              %%Keyword: aRV <basic>
              <HELP>
              Read the real number R value around which the radial "variable" 
              used for numerically solving the differential equation, is
              centered.
              </HELP>
              </KEYWORD>

:kword:`EPS`
  The real number :math:`\epsilon` value indicating the convergence tolerance when numerically solving the differential equation.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="EPS" KIND="REAL" LEVEL="BASIC">
              %%Keyword: EPS <basic>
              <HELP>
              Read the real number epsilon value indicating the convergence 
              tolerance when numerically solving the differential equation.
              </HELP>
              </KEYWORD>

:kword:`NTP`
  The integer indicating the number of turning points.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="NTP" KIND="INT" LEVEL="BASIC">
              %%Keyword: NTP <basic>
              <HELP>
              Read the integer indicating the number of turning points when 
              providing a pointwise potential in the input file.
              </HELP>
              </KEYWORD>

:kword:`LPPOt`
  The integer indicating how often to print the potential and its first two derivatives (they will all be printed to Channel 6 at every (LPPOT)th point if LPPOT > 0, and only the potential will be printed in condensed format to Channel 8 at every \|LPPOT\|th point if LPPOT < 0). 

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="LPPOT" KIND="INT" LEVEL="UNDOCUMENTED">
              %%Keyword: LPPOt <undocumented>
              <HELP>
              Read the integer indicating how often to print the potential 
              and its first two derivatives.
              </HELP>
              </KEYWORD>

:kword:`IOMEg`
  The integer angular momentum quantum number :math:`\Omega`.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME=IOMEG" KIND="INT" LEVEL="UNDOCUMENTED">
              %%Keyword: IOMEg <undocumented>
              <HELP>
              Read the integer angular momentum quantum number Omega.
              </HELP>
              </KEYWORD>

:kword:`VLIM`
  The real number indicating the limit of the potential :math:`V(R)` as :math:`R\rightarrow \infty`.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="VLIM" KIND="REAL" LEVEL="BASIC">
              %%Keyword: VLIM <basic>
              <HELP>
              Read the real number indicating the limit of the potential 
              V(R) as R -> infinity.
              </HELP>
              </KEYWORD>

:kword:`IPOTl`
  The integer indicating the form of the analytic potential. Choose IPOTL = 1 for a Lennard-Jones potential, IPOTL = 3 for an EMO (extended Morse oscillator), IPOTL = 4 for the MLR (Morse/Long-range) potential.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="IPOTL" KIND="INT" LEVEL="UNDOCUMENTED">
              %%Keyword: IPOTl <undocumented>
              <HELP>
              Read the integer indicating the form of the analytic potential
              being used.
              </HELP>
              </KEYWORD>

:kword:`PPAR`
  The integer power :math:`p` used in an MLR potential. 

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="PPAR" KIND="INT" LEVEL="BASIC">
              %%Keyword: PPAR <basic>
              <HELP>
              Read the integer power p used in an MLR potential
              </HELP>
              </KEYWORD>

:kword:`QPAR`
  The integer power :math:`q` used in an MLR potential. 

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="QPAR" KIND="INT" LEVEL="BASIC">
              %%Keyword: QPAR <basic>
              <HELP>
              Read the integer power q used in an MLR potential
              </HELP>
              </KEYWORD>

:kword:`NSR`
  The integer order of the polynomial function in an MLR potential's exponent, for the short-range (SR) part of the potential. 

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="NSR" KIND="INT" LEVEL="UNDOCUMENTED">
              %%Keyword: NSR <undocumented>
              <HELP>
              Read the integer order of the polynomial function in an MLR
              potential's exponent, for the short-range (SR) part of the potential. 
              </HELP>
              </KEYWORD>

:kword:`NLR`
  The integer order of the polynomial function in an MLR potential's exponent, for the long-range (LR) part of the potential. 

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="NLR" KIND="INT" LEVEL="UNDOCUMENTED">
              %%Keyword: NLR <undocumented>
              <HELP>
              Read the integer order of the polynomial function in an MLR 
              potential's exponent, for the long-range (LR) part of the potential. 
              </HELP>
              </KEYWORD>

:kword:`IBOB`
  The integer flag specifying whether or not to include (IBOB>0) or exclude (IBOB <= 0) Born-Oppenheimer Breakdown functions.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="IBOB" KIND="INT" LEVEL="UNDOCUMENTED">
              %%Keyword: IBOB <undocumented>
              <HELP>
              Read the integer flag specifying whether or not to include (IBOB>0)
              or exclude (IBOB <= 0) Born-Oppenheimer Breakdown functions.
              </HELP>
              </KEYWORD>

:kword:`DSCM`
  The real number indicating the :math:`\mathfrak{D}_e` value (the "depth at equilibrium" for the potential).

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="DSCM" KIND="REAL" LEVEL="BASIC">
              %%Keyword: DSCM <basic>
              <HELP>
              Read the real number indicating the De value (the "depth at
              equilibrium" for the potential).
              </HELP>
              </KEYWORD>

:kword:`REQ`
  The real number indicating the :math:`R_e` value (the equilibrium internuclear distance)

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="REQ" KIND="REAL" LEVEL="BASIC">
              %%Keyword: REQ <basic>
              <HELP>
              The real number indicating the R_e value (the equilibrium internuclear
              distance).
              </HELP>
              </KEYWORD>

:kword:`RREF`
  The reasl number indicating the :math:`R_ref` value (the reference distance for the MLR model).
  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="RREF" KIND="REAL" LEVEL="BASIC">
              %%Keyword: RREF <basic>
              <HELP>
              Read the real number indicating the "reference distance"a round which
              the MLR model is "centered".
              </HELP>
              </KEYWORD>

:kword:`NCMM`
  Integer indicating the number of long-range terms used in the MLR model (e.g. if using C6,C8,C10, then NCMM=3).

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="NCMM" KIND="INT" LEVEL="UNDOCUMENTED">
              %%Keyword: NCMM <undocumented>
              <HELP>
              Read the integer indicating how many long-range terms to include
              in the MLR potential.
              </HELP>
              </KEYWORD>

:kword:`IVSR`
  Integer indicating the power used in the damping function for the MLR model.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="IVSR" KIND="INT" LEVEL="UNDOCUMENTED">
              %%Keyword: IVSR  <undocumented>
              <HELP>
              Read the integer indicating how often to print the potential 
              and its first two derivatives.
              </HELP>
              </KEYWORD>

:kword:`IDSTt`
  Integer indicating the type of damping function used. Choose "1" for the Douketis-Scoles-type (DS) function, and "2" for the Tang-Toonies (TS) function.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="IDSTT" KIND="INT" LEVEL="UNDOCUMENTED">
              %%Keyword: IDSTt <undocumented>
              <HELP>
              Read the integer indicating which type of damping function to use in
              the MLR model.
              </HELP>
              </KEYWORD>

:kword:`RHOAb`
  Real number indicating the :math:`\rho_{AB}` parameter for the damping function in an MLR model.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="RHOAB" KIND="REAL" LEVEL="UNDOCUMENTED">
              %%Keyword: RHOAb <undocumented>
              <HELP>
              Read the real number indicating the value of the rho_AB damping function 
              parameter for an MLR model.
              </HELP>
              </KEYWORD>

:kword:`MMLR`
  Integer array containing NCMM elements, which indicate the inverse powers of the long-range terms in the MLR model. For example, if using C6,C8,C10, then MMLR = 6 8 10. If using C4,C6,C8 (for example, for the potential between a neutral atom and an ion) then use MMLR = 4 6 8.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="MMLR" KIND="INTS" SIZE="3" LEVEL="UNDOCUMENTED">
              %%Keyword: MMLR <undocumented>
              <HELP>
              Read the integer array indicating the values of the inverse-powers for 
              the long-range tail of an MLR model.
              </HELP>
              </KEYWORD>

:kword:`CMM`
  Real number array containing NCMM elements, which indicate the coefficients of the inverse powers of the long-range terms in the MLR model. For example, if using C6,C8
,C10, then CMM = C6 C8 C10. If using C4,C6,C8 (for example, for the potential between a neutral atom and an ion) then use MMLR = C4 C6 C8.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="CMM" KIND="REALS" SIZE="3" LEVEL="UNDOCUMENTED">
              %%Keyword: CMM <undocumented>
              <HELP>
              Read the real-number array indicating the values of the coefficients of
              the inverse-powers for the long-range tail of an MLR model.
              </HELP>
              </KEYWORD>

:kword:`PARM`
  Real number array containing NLR elements, which indicate the exponent expansion coefficients for the MLR model.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="PARM" KIND="REALS" SIZE="4" LEVEL="BASIC">
              %%Keyword: PARM <BASIC>
              <HELP>
              Read the real-number array indicating the values of the exponent
              expansion coefficients for the MLR model.
              </HELP>
              </KEYWORD>

:kword:`NLEV1`
  Integer indicating the number of rovibrational levels to seek. If negative, the program will try to automatically find all levels from :math:`v=0` to `v=-|\textrm{NLEV1}|`. 

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="NLEV1" KIND="INT" LEVEL="BASIC">
              %%Keyword: NLEV1 <BASIC>
              <HELP>
              Read the integer indicating the number of rovibrational levels
              to find.
              </HELP>
              </KEYWORD>

:kword:`AUTO1`
  Integer indicating whether or not to automatically generate trial energies for each vibrational level. If > 0, the trial energies are generated, wheras if <= 0, then the user can provide trial energies manually.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="AUTO1" KIND="INT" LEVEL="UNDOCUMENTED">
              %%Keyword: AUTO1 <undocumented>
              <HELP>
              Read the integer indicating whether or not to automatically 
              generate trial energies for each vibrational level sought.
              </HELP>
              </KEYWORD>

:kword:`LCDC`
  Integer indicating whether or not to calculate inertial rotational constants: :math:`B_v`, and the first six centrifugal distortion constants: :math:`-D_v,H_v,L_v,M_v,N_v,O_v`. If >0, then these are calculated, and otherwise they are not. 

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="LCDC" KIND="INT" LEVEL="BASIC">
              %%Keyword: LCDC <basic>
              <HELP>
              Integer indicating whether or not to calculate Bv,-Dv,Hv,Lv,Mv,Nv,Ov.
              </HELP>
              </KEYWORD>

:kword:`LXPCt`
  Integer indicating whether or not to calculate expectation values or matrix elements using the ro-vibrational wavefunctions obtained from solving the Schroedinger equation. If =0, no expectation values or matrix elements are calculated, and otherwise they are.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="LXPCT" KIND="INT" LEVEL="UNDOCUMENTED">
              %%Keyword: LXPCt <undocumented>
              <HELP>
              Read the integer indicating whether or not to print expectation values
              or matrix elements.
              </HELP>
              </KEYWORD>

:kword:`NJM`
  Integer indicating how many rotational levels (and expectation values, if LXPCT>0) to find for each vibrational level found. 

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="NJM" KIND="INT" LEVEL="UNDOCUMENTED">
              %%Keyword: NJM <undocumented>
              <HELP>
              Read the integer indicating how many rotational levels to find for 
              each vibrational level found.
              </HELP>
              </KEYWORD>


:kword:`LPRWf`
  Integer indicating whether or not to print the ro-vibrational wavefunction levels at every LPRWF'th mesh point. If =0, no wavefunction is printed.

  .. xmldoc:: <KEYWORD MODULE="LEVEL" NAME="LPRWF" KIND="INT" LEVEL="UNDOCUMENTED">
              %%Keyword: LPRWf <undocumented>
              <HELP>
              Read the ineger indicating whether or not to print the wavefunction.
              </HELP>
              </KEYWORD>




Input example
.............

::

  &LEVEL
    IAN1 = 3
    IMN1 = 6
    IAN2 = 3
    IMN2 = 6
    CHARGE = 0
    NUMPOT = 1
    RH = 0.0005
    RMIN = 0.125
    PRV = 1
    ARV = 5.0d0
    EPS = 2.d-10
    NTP = -1
    LPPOT = 0
    IOMEG1 = 0
    VLIM = 0.0d0
    IPOTL = 4
    PPAR = 5
    QPAR = 3
    NSR = 3
    NLR = 3
    IBOB = -1
    DSCM = 3.337678701485D+02
    REQ = 4.170010583477D+00
    RREF = 8.0d0
    NCMM = 3
    IVSR = -2
    TDSTT = 1
    RHOAB = 0.54d0
    MMLR = 6 8 10
    CMM = 6.719000000d+06 1.126350000d+08  2.786940000d+09
    PARM = -5.156803528943D-01 -9.585070416286D-02 1.170797201140D-01 -2.282814434665D-02
    NLEV1 = -999
    AUTO1 = 1
    LCDC = 2
    LXPCT = 0
    NJM = 0
    JDJR = 1
    LPRWF = 0

**Comments**: The vibrational-rotation spectrum for the :math:`1^3\Sigma_u(a)` state of
 :math:`^{(6,6)}\ce{Li2}` will be computed using the MLR potential given in the input. 

.. xmldoc:: </MODULE>
