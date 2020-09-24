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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine CoW(Coor,CoF,W,nAtom,T)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January '91                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Coor(3,nAtom), CoF(3), T, W(nAtom)
*
      iRout = 140
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.99) Then
         Call RecPrt(' In CoW: Coor',' ',Coor,3,nAtom)
         Call RecPrt(' In CoW: W',' ',W,nAtom,1)
      End If
      T = Zero
      Do 10 iAtom = 1, nAtom
         T=T+W(iAtom)
 10   Continue
      Do 20 iCar = 1, 3
         CoF(iCar) = Zero
         Do 21 iAtom = 1, nAtom
            CoF(iCar)=CoF(iCar) + Coor(iCar,iAtom)*W(iAtom)
 21      Continue
         If (T.ne.Zero) Then
            CoF(iCar)= CoF(iCar) / T
         Else
            CoF(iCar)= Zero
         End If
 20   Continue
      If (iPrint.ge.99) Then
         Call RecPrt(' In CoW: CoF',' ',CoF,1,3)
         Call RecPrt(' In CoW: T',' ',[T],1,1)
      End If
*
      Return
      End
