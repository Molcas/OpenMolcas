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
      SubRoutine ModU2(U2,mT,nRys,ZEInv)
************************************************************************
*                                                                      *
* Object: precompute u2/(zeta+eta)                                     *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             May '90                                                  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Real*8 U2(nRys,mT), ZEInv(mT)
*
      iRout = 255
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.99) Then
         Call RecPrt(' In ModU2: U2',' ',U2,nRys,mT)
         Call RecPrt(' In ModU2: ZEInv',' ',ZEInv,1,mT)
      End If
*
      If (nRys.gt.1) Then
         Do 11 iT = 1, mT
            Do 21 iRys = 1, nRys
                  U2(iRys,iT) = U2(iRys,iT) * ZEInv(iT)
 21         Continue
 11      Continue
      Else
         Do 31 iT = 1, mT
            U2(1,iT) = U2(1,iT) * ZEInv(iT)
 31      Continue
      End If
*
      Return
      End
