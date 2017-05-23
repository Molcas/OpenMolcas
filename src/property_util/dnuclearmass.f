************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2000, Per-Olof Widmark                                 *
************************************************************************
      Real*8 Function dNuclearMass(Z,A)
************************************************************************
*                                                                      *
* Routine: dNuclearMass                                                *
* Purpose: This is a wrapper for dxNuclearMass to compute the mass for *
*          isotope with Z protons and A-Z neutrons.                    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden.                                    *
* Written: March 2000                                                  *
* History: none.                                                       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Parameters:                                                          *
* Z  - The nuclear charge for the isotope.                             *
* A  - The number of nucleons for the isotope.                         *
*                                                                      *
************************************************************************
      Implicit  None
*----------------------------------------------------------------------*
* Dummy parameters.                                                    *
*----------------------------------------------------------------------*
      Integer   Z
      Integer   A
*----------------------------------------------------------------------*
* Local variables.                                                     *
*----------------------------------------------------------------------*
      Real*8    Mass
      Integer   Rc
      Integer   Opt
*----------------------------------------------------------------------*
* External references.                                                 *
*----------------------------------------------------------------------*
      Real*8    dxNuclearMass
      External  dxNuclearMass
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Rc=0
      Opt=0
      Mass=dxNuclearMass(Z,A,Rc,Opt)
      If(Rc.ne.0) Then
         Call SysAbendMsg('dnuclearmass','fail to get mass',' ')
      End If
*----------------------------------------------------------------------*
* Done.                                                                *
*----------------------------------------------------------------------*
      dNuclearMass=Mass
      Return
      End
