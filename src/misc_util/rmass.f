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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      Function rMass(nAtom)
************************************************************************
*                                                                      *
* Object: to return the mass of the nucleus as a function of the       *
*         atomic number, nAtom. The mass is that one of the most       *
*         abundant isotope. In the case there is not stable isotope    *
*         we select the one with the longest lifetime.                 *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
************************************************************************
      Use Isotopes
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "constants2.fh"
      Real*8 rMass
*
      rMass=Zero
      If (nAtom.gt.MaxAtomNum) Then
*        Write (6,*) ' Weight for this atom is not listed!'
*        Write (6,*) ' Mass set to 2.6 times atom number'
         rMass = 2.6D0 * DBLE(nAtom) * uToau
      Else If (nAtom.eq.0) Then
*        Write (6,*) ' Weight for this atom is meaningless!'
*        Write (6,*) ' Mass set to 0.0'
      Else If (nAtom.lt.0) Then
         rMass=1.0D99 * uToau
      Else
         isnx=0
         Call Isotope(isnx,nAtom,rMass)
      End If
      Return
      End
