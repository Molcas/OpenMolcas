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
      Subroutine sssp(EFInt,Zeta,nZeta,P,lP,rKappAB,A,B,
     &                       Eta, nEta,Q,lQ,rKappCD,C,D,
     &                CoorAC,TMax,
     &                iPntr,nPntr,x0,nMax,W6,W5,W4,W3,W2,W1,W0,
     &                                    R6,R5,R4,R3,R2,R1,R0,
     &                ddx,HerW,HerR2,IsChi,ChiI2)
************************************************************************
*                                                                      *
* Object: to compute the primitive integrals of type (ss|sp).          *
*                                                                      *
*  Author:    Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN. 1994                                    *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 EFInt(nZeta,nEta,3), Zeta(nZeta), Eta(nEta),
     &       CoorAC(3,2),
     &       P(lP,3), Q(lQ,3), A(3), B(3), C(3), D(3),
     &       rKappAB(nZeta), rKappCD(nEta),
     &       x0(nMax), W6(nMax), W5(nMax),
     &       W4(nMax), W3(nMax), W2(nMax), W1(nMax), W0(nMax),
     &       R6(nMax), R5(nMax),
     &       R4(nMax), R3(nMax), R2(nMax), R1(nMax), R0(nMax)
      Integer iPntr(nPntr)
      Logical ABeqCD, EQ
*
*
      xdInv=One/ddx
      dddx = ddx/10d0 + ddx
*
      ABeqCD = EQ(A,B) .and. EQ(A,C) .and. EQ(A,D)
      If (ABeqCD) Go To 300
      If (EQ(C,D)) Go To 200
*
*-----ABCD case
*
      Do 10 iEta = 1, nEta
         Do 20 iZeta = 1, nZeta
            ZEInv = One/(Eta(iEta)+Zeta(iZeta)
     >              +(Eta(iEta)*Zeta(iZeta)*ChiI2)*Dble(IsChi))
            ZE  = Zeta(iZeta)*Eta(iEta)
            rho = ZE*ZEInv
            PQx = P(iZeta,1)-Q(iEta,1)
            PQy = P(iZeta,2)-Q(iEta,2)
            PQz = P(iZeta,3)-Q(iEta,3)
            PQ2 = PQx**2 + PQy**2 + PQz**2
            T = rho*PQ2
            If (T.lt.TMax) Then
               n = iPntr(Int((T+dddx)*xdInv))
               z = T - x0(n)
               w =(((((W6(n)*z+W5(n))*z+W4(n))*z+W3(n))*z+W2(n))
     &           *z+W1(n))*z+w0(n)
               r =(((((R6(n)*z+R5(n))*z+R4(n))*z+R3(n))*z+R2(n))
     &           *z+R1(n))*z+R0(n)
               Zu2 = r * (Zeta(iZeta)*ZEInv)
               PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv) * w
            Else
               Zu2 = HerR2 / (Eta(iEta)*PQ2)
               PreFct = rKappCD(iEta) * rKappAB(iZeta)
     &                * HerW / Sqrt(ZE*PQ2)
            End If
            QCPQx = Q(iEta,1) - CoorAC(1,2) + Zu2 * PQx
            QCPQy = Q(iEta,2) - CoorAC(2,2) + Zu2 * PQy
            QCPQz = Q(iEta,3) - CoorAC(3,2) + Zu2 * PQz
            EFInt(iZeta,iEta,1) = PreFct * QCPQx
            EFInt(iZeta,iEta,2) = PreFct * QCPQy
            EFInt(iZeta,iEta,3) = PreFct * QCPQz
 20      Continue
 10   Continue
      Go To 99
*
*-----ABCC case
*
 200  Continue
      Do 11 iEta = 1, nEta
         Do 21 iZeta = 1, nZeta
            ZEInv = One/(Eta(iEta)+Zeta(iZeta)
     >              +(Eta(iEta)*Zeta(iZeta)*ChiI2)*Dble(IsChi))
            ZE  = Zeta(iZeta)*Eta(iEta)
            rho = ZE*ZEInv
            PQx =  (P(iZeta,1)-CoorAC(1,2))
            PQy =  (P(iZeta,2)-CoorAC(2,2))
            PQz =  (P(iZeta,3)-CoorAC(3,2))
            PQ2 = PQx**2 + PQy**2 + PQz**2
            T = rho*PQ2
            If (T.lt.TMax) Then
               n = iPntr(Int((T+dddx)*xdInv))
               z = T - x0(n)
               w =(((((W6(n)*z+W5(n))*z+W4(n))*z+W3(n))*z+W2(n))
     &           *z+W1(n))*z+w0(n)
               r =(((((R6(n)*z+R5(n))*z+R4(n))*z+R3(n))*z+R2(n))
     &           *z+R1(n))*z+R0(n)
               Zu2 = r * (Zeta(iZeta)*ZEInv)
               PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv) * w
            Else
               Zu2 = HerR2 / (Eta(iEta)*PQ2)
               PreFct = rKappCD(iEta) * rKappAB(iZeta)
     &                * HerW / Sqrt(ZE*PQ2)
            End If
            EFInt(iZeta,iEta,1) = PreFct * Zu2 * PQx
            EFInt(iZeta,iEta,2) = PreFct * Zu2 * PQy
            EFInt(iZeta,iEta,3) = PreFct * Zu2 * PQz
 21      Continue
 11   Continue
      Go To 99
*
*-----CCCC case
*
 300  Continue
      Do iEta = 1, nEta
         Do iZeta = 1, nZeta
            EFInt(iZeta,iEta,1) = Zero
            EFInt(iZeta,iEta,2) = Zero
            EFInt(iZeta,iEta,3) = Zero
         End Do
      End Do
*
 99   Continue
*
      Return
      End
