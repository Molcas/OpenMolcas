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
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      SubRoutine PrePre_g(nZeta,nEta,mZeta,mEta,lZeta,lEta,
     &                  Data1,Data2,PreScr,CutGrd)
************************************************************************
*                                                                      *
* Object: to preprescreen the integral derivatives.                    *
*                                                                      *
*   nZeta, nEta : unpartioned length of primitives.                    *
*                                                                      *
*   mZeta, mEta : section length due to partioning. These are usually  *
*                 equal to nZeta and nEta.                             *
*                                                                      *
* Called from: Twoel                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             July '92.                                                *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8
     &       Data1(nZeta,8), Data2(nEta,8)
      Logical PreScr
#include "print.fh"
#include "real.fh"
*
      iRout = 180
      iPrint = nPrint(iRout)
*     iQ = 0
*     Call qEnter('PrePre')
      If (iPrint.ge.99) Then
         Call RecPrt(' Data1',' ',Data1,nZeta,8)
         Call RecPrt(' Data2',' ',Data2,nEta ,8)
      End If
*
      tOne = One
      ip=1
*
*-----Preprescanning
*
      lZeta=mZeta
      lEta =mEta
      rKabMx=Zero
      ZetaMx=Zero
      rKabMn=1.0D+72
      ZetaMn=Zero
      Do 700 iZeta = 1, mZeta
         If (Data1(iZeta,2).gt.rKabMx) Then
            rKabMx=Data1(iZeta,2)
            ZetaMx=Data1(iZeta,1)
         End If
         If (Data1(iZeta,2).lt.rKabMn) Then
            rKabMn=Data1(iZeta,2)
            ZetaMn=Data1(iZeta,1)
         End If
 700  Continue
      rKcdMx=Zero
      EtaMx =Zero
      rKcdMn=1.0D+72
      EtaMn =Zero
      Do 701 iEta = 1, mEta
         If (Data2(iEta,2).gt.rKcdMx) Then
            rKcdMx=Data2(iEta,2)
            EtaMx =Data2(iEta,1)
         End If
         If (Data2(iEta,2).lt.rKcdMn) Then
            rKcdMn=Data2(iEta,2)
            EtaMn =Data2(iEta,1)
         End If
 701  Continue
      PreScr=.True.
      PreMax=rKabMx*rKcdMx*Sqrt(tOne/(ZetaMx+EtaMx))
      PreMin=rKabMn*rKcdMn*Sqrt(tOne/(ZetaMn+EtaMn))
      If (PreMin.gt.CutGrd) PreScr=.False.
      If (PreMax.lt.1.0D-4*CutGrd) Then
         lZeta = 0
         lEta = 0
      End If
*
*     Call qExit('PrePre')
      Return
      End
