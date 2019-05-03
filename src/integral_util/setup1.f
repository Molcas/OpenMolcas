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
      SubRoutine Setup1(ExpA,nPrim,ExpB,mPrim,A,B,rKappa,Pcoor,ZInv)
************************************************************************
*                                                                      *
*     Object : to compute some data which is needed for the one-       *
*              electron integrals.                                     *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : RecPrt                                                  *
*              GetMem                                                  *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             January '90                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 ExpA(nPrim), ExpB(mPrim), rKappa(nPrim,mPrim),
     &       Pcoor(nPrim,mPrim,3),  ZInv(nPrim,mPrim), A(3), B(3)
*
      iRout = 114
      iPrint = nPrint(iRout)
*     Call qEnter('SetUp1')
      ab  = (A(1)-B(1))**2 + (A(2)-B(2))**2 + (A(3)-B(3))**2
      If (ab.ne.Zero) Then
      Do 10 iPrim = 1, nPrim
         Do 20 jPrim = 1, mPrim
            rKappa(iPrim,jPrim) = Exp(- ExpA(iPrim) * ExpB(jPrim) * ab *
     &                            ZInv(iPrim,jPrim))
            Pcoor(iPrim,jPrim,1)=(ExpA(iPrim)*A(1)+ExpB(jPrim)*B(1)) *
     &                            ZInv(iPrim,jPrim)
            Pcoor(iPrim,jPrim,2)=(ExpA(iPrim)*A(2)+ExpB(jPrim)*B(2)) *
     &                            ZInv(iPrim,jPrim)
            Pcoor(iPrim,jPrim,3)=(ExpA(iPrim)*A(3)+ExpB(jPrim)*B(3)) *
     &                            ZInv(iPrim,jPrim)
 20      Continue
 10   Continue
      Else
        call dcopy_(nPrim*mPrim,[One],0,rKappa,1)
        call dcopy_(nPrim*mPrim,A(1),0,Pcoor(1,1,1),1)
        call dcopy_(nPrim*mPrim,A(2),0,Pcoor(1,1,2),1)
        call dcopy_(nPrim*mPrim,A(3),0,Pcoor(1,1,3),1)
      End If
      If (iPrint.ge.99) Then
         Call RecPrt(' *** Kappa ***',' ',rKappa, nPrim, mPrim)
         Call RecPrt(' ***   Px  ***',' ',Pcoor(1,1,1),nPrim,mPrim)
         Call RecPrt(' ***   Py  ***',' ',Pcoor(1,1,2),nPrim,mPrim)
         Call RecPrt(' ***   Pz  ***',' ',Pcoor(1,1,3),nPrim,mPrim)
      End If
*     Call GetMem('SetUp1','CHECK','REAL',iDum,iDum)
*
*     Call qExit('SetUp1')
      Return
      End
