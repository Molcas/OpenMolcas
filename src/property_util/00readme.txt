********************************************************************************
*                                                                              *
*                              Property utilities                              *
*                                                                              *
********************************************************************************

This directory contains utility routines for properties of various kinds. The
routines are named in a style resembling the linpack conventions: the first
letter denote the data type the routine process:

o i -- integer data type
o s -- real*4 data type
o d -- real*8 data type
o c -- complex*4 data type
o z -- complex*8 data type

The routines come in in two versions, one with a minimum of parameters and a
default behavior, one with return code and options added to the parameter list.
The latter is denoted with the letter x after the first type character, for
example dNuclearMass and dxNuclearMass do the same work but the latter gives
more control to the user.

--------------------------------------------------------------------------------

Real*8 Function dNuclearMass(Z,A)
Real*8 Function dxNuclearMass(Z,A,Rc,Opt)

This routine returns the nuclear mass of a given isotope number.

Parameters:
Z   - Integer, input. The charge of the nucleus, i.e. the number of protons.
A   - Integer, input. The mass number of the nucleus, i.e. the number of protons
      plus the number of neutrons.
Rc  - Integer, output. The return code. This is zero if successful and nonzero if
      not.
Opt - Integer, input. Options switch to modify the default behavior. Option 0
      (zero) gives default behavior.

This routine returns the nuclear mass of a given isotope number. The mass
is returned in atomic units (electron masses), not atomic mass units (1/12*C12
masses). Most isotopes have their masses tabulated, but a few exotic are not
due to lack of experimental data or other reasons. For such isotopes, the
semi-empirical mass formula is used and a warning is printed.

--------------------------------------------------------------------------------

Integer Function iMostAbundantIsotope(Z)
Integer Function ixMostAbundantIsotope(Z,Rc,Opt)

This routine returns the mass number of the most abundant isotope for a given
nuclear charge.

Parameters:
Z   - Integer, input. The charge of the nucleus, i.e. the number of protons.
Rc  - Integer, output. The return code. This is zero if successful and nonzero if
      not.
Opt - Integer, input. Options switch to modify the default behavior. Option 0
      (zero) gives default behavior.

This routine returns the mass number (number of protons plus neutrons) of the
most abundant isotope for a given nuclear charge. In the case of radioactive
isotopes with no natural occurance the most stable isotope is returned.

--------------------------------------------------------------------------------

Integer Function iNuclearChargeFromSymbol(Symbol)
Integer Function ixNuclearChargeFromSymbol(Symbol,Rc,Opt)

This routine returns the nuclear charge of an atom corresponding to the
chemical symbol.

Parameters:
Symbol - Character*(*), input. The chemical symbol of the atom for which you
         want the nuclear charge, 'H', 'He', 'Li', ...
Rc     - Integer, output. The return code. This is zero if successful and nonzero
         if not.
Opt    - Integer, input. Options switch to modify the default behavior. Option
         0 (zero) gives default behavior.

This routine returns the nuclear charge of an atom corresponding to the
chemical symbol. The case of the chemical symbol does not matter, 'He', 'he',
'HE' and 'hE' all return charge 2. All symbols up to 'Og' (118) are tabulated.

--------------------------------------------------------------------------------

