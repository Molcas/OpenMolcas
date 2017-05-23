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
* Copyright (C) 1994, Roland Lindh                                     *
************************************************************************
      Subroutine ssss(EFInt,Zeta,nZeta,P,lP,rKappAB,A,B,
     &                       Eta, nEta,Q,lQ,rKappCD,C,D,TMax,
     &                iPntr,nPntr,x0,nMax,W6,W5,W4,W3,W2,W1,W0,ddx,HerW,
     &                IsChi,ChiI2)
************************************************************************
*                                                                      *
* Object: to compute the primitive integrals of type (ss|ss).          *
*                                                                      *
* Called from: vRys                                                    *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*                                                                      *
*  Author:    Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN. 1994                                    *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 EFInt(nZeta,nEta), Zeta(nZeta), Eta(nEta),
     &       P(lP,3), Q(lQ,3), A(3), B(3), C(3), D(3),
     &       rKappAB(nZeta), rKappCD(nEta),
     &       x0(nMax), W6(nMax), W5(nMax),
     &       W4(nMax), W3(nMax), W2(nMax), W1(nMax), W0(nMax)
      Integer iPntr(nPntr)
      Logical ABeqCD, EQ
*
*     Call qEnter('ssss')
*
      xdInv=One/ddx
      dddx = ddx/10d0 + ddx
*
      ABeqCD = EQ(A,B) .and. EQ(A,C) .and. EQ(A,D)
      If (ABeqCD) Go To 100
*
      Do 10 iEta = 1, nEta
         Do 20 iZeta = 1, nZeta
            ZEInv = One/(Eta(iEta)+Zeta(iZeta)
     >              +(Eta(iEta)*Zeta(iZeta)*ChiI2)*Dble(IsChi))
            ZE  = Zeta(iZeta)*Eta(iEta)
            rho = ZE*ZEInv
            PQ2 = (P(iZeta,1)-Q(iEta,1))**2
     &          + (P(iZeta,2)-Q(iEta,2))**2
     &          + (P(iZeta,3)-Q(iEta,3))**2
            T = rho*PQ2
            If (T.lt.TMax) Then
               n = iPntr(Int((T+dddx)*xdInv))
               z = T - x0(n)
               w =(((((W6(n)*z+W5(n))*z+W4(n))*z+W3(n))*z+W2(n))
     &           *z+W1(n))*z+w0(n)
               EFInt(iZeta,iEta) = rKappCD(iEta) * rKappAB(iZeta)
     &                           * Sqrt(ZEInv) * w
            Else
               EFInt(iZeta,iEta) = rKappCD(iEta) * rKappAB(iZeta)
     &                           * HerW * Sqrt(One/(ZE*PQ2))
            End If
 20      Continue
 10   Continue
      Go To 99
*
 100  Continue
      z = - x0(1)
      w =(((((W6(1)*z+W5(1))*z+W4(1))*z+W3(1))*z+W2(1))
     &   *z+W1(1))*z+w0(1)
      Do 11 iEta = 1, nEta
         Do 21 iZeta = 1, nZeta
            ZEInv = One/(Eta(iEta)+Zeta(iZeta)
     >              +(Eta(iEta)*Zeta(iZeta)*ChiI2)*Dble(IsChi))
            EFInt(iZeta,iEta) = rKappCD(iEta) * rKappAB(iZeta)
     &                        * Sqrt(ZEInv) * w
 21      Continue
 11   Continue
*
 99   Continue
*
*     Call qExit('ssss')
      Return
      End
