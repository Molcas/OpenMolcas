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
      Subroutine psps(EFInt,Zeta,nZeta,P,lP,rKappAB,A,B,
     &                       Eta, nEta,Q,lQ,rKappCD,C,D,
     &                CoorAC,TMax,
     &                iPntr,nPntr,x0,nMax,CW6,CW5,CW4,CW3,CW2,CW1,CW0,
     &                                    CR6,CR5,CR4,CR3,CR2,CR1,CR0,
     &                ddx,HerW,HerR2,IsChi,ChiI2)
************************************************************************
*                                                                      *
* Object: to compute the primitive integrals of type (ps|ps).          *
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
      Real*8 EFInt(nZeta,nEta,3,3), Zeta(nZeta), Eta(nEta),
     &       CoorAC(3,2),
     &       P(lP,3), Q(lQ,3), A(3), B(3), C(3), D(3),
     &       rKappAB(nZeta), rKappCD(nEta),
     &       x0(nMax),
     &       CW6(nMax,2), CW5(nMax,2), CW4(nMax,2), CW3(nMax,2),
     &       CW2(nMax,2), CW1(nMax,2), CW0(nMax,2),
     &       CR6(nMax,2), CR5(nMax,2), CR4(nMax,2), CR3(nMax,2),
     &       CR2(nMax,2), CR1(nMax,2), CR0(nMax,2),
     &       HerW(2), HerR2(2)
      Integer iPntr(nPntr)
      Logical ABeqCD, EQ
*
*     Call qEnter('psps')
*
      xdInv=One/ddx
      dddx = ddx/10d0 + ddx
*
      ABeqCD = EQ(A,B) .and. EQ(A,C) .and. EQ(A,D)
      If (           ABeqCD            ) Go To 100
      If (     EQ(A,B).and..Not.EQ(C,D)) Go To 200
      If (.Not.EQ(A,B).and.     EQ(C,D)) Go To 300
      If (     EQ(A,B).and.     EQ(C,D)) Go To 400
*
*-----ABCD case
*
      Do 10 iEta = 1, nEta
         Do 20 iZeta = 1, nZeta
            ZEInv = One/(Eta(iEta)+Zeta(iZeta)
     >              +(Eta(iEta)*Zeta(iZeta)*ChiI2)*Dble(IsChi))
            rho = Zeta(iZeta)*Eta(iEta)*ZEInv
            PQx = P(iZeta,1)-Q(iEta,1)
            PQy = P(iZeta,2)-Q(iEta,2)
            PQz = P(iZeta,3)-Q(iEta,3)
            T = rho * (PQx**2 + PQy**2 + PQz**2)
            If (T.lt.TMax) Then
               n = iPntr(Int((T+dddx)*xdInv))
               z = T - x0(n)
               w1=(((((CW6(n,1)*z+CW5(n,1))*z+CW4(n,1))*z+CW3(n,1))*z+
     &            CW2(n,1))*z+CW1(n,1))*z+Cw0(n,1)
               w2=(((((CW6(n,2)*z+CW5(n,2))*z+CW4(n,2))*z+CW3(n,2))*z+
     &            CW2(n,2))*z+CW1(n,2))*z+Cw0(n,2)
               r1=(((((CR6(n,1)*z+CR5(n,1))*z+CR4(n,1))*z+CR3(n,1))*z+
     &            CR2(n,1))*z+CR1(n,1))*z+CR0(n,1)
               r2=(((((CR6(n,2)*z+CR5(n,2))*z+CR4(n,2))*z+CR3(n,2))*z+
     &            CR2(n,2))*z+CR1(n,2))*z+CR0(n,2)
            Else
               ai = 1.0D0/T
               si = Sqrt(ai)
               w1= HerW(1)*si
               w2= HerW(2)*si
               r1= HerR2(1)*ai
               r2= HerR2(2)*ai
            End If
            Eu21 =   Eta(iEta) * (r1*ZEInv)
            Eu22 =   Eta(iEta) * (r2*ZEInv)
            Zu21 = Zeta(iZeta) * (r1*ZEInv)
            Zu22 = Zeta(iZeta) * (r2*ZEInv)
            PAQPx1 = (P(iZeta,1) - CoorAC(1,1)) - Eu21 * PQx
            PAQPx2 = (P(iZeta,1) - CoorAC(1,1)) - Eu22 * PQx
            PAQPy1 = (P(iZeta,2) - CoorAC(2,1)) - Eu21 * PQy
            PAQPy2 = (P(iZeta,2) - CoorAC(2,1)) - Eu22 * PQy
            PAQPz1 = (P(iZeta,3) - CoorAC(3,1)) - Eu21 * PQz
            PAQPz2 = (P(iZeta,3) - CoorAC(3,1)) - Eu22 * PQz
            QCPQx1 = (Q(iEta,1)  - CoorAC(1,2)) + Zu21 * PQx
            QCPQx2 = (Q(iEta,1)  - CoorAC(1,2)) + Zu22 * PQx
            QCPQy1 = (Q(iEta,2)  - CoorAC(2,2)) + Zu21 * PQy
            QCPQy2 = (Q(iEta,2)  - CoorAC(2,2)) + Zu22 * PQy
            QCPQz1 = (Q(iEta,3)  - CoorAC(3,2)) + Zu21 * PQz
            QCPQz2 = (Q(iEta,3)  - CoorAC(3,2)) + Zu22 * PQz
            B001 = Half * (r1*ZEInv)
            B002 = Half * (r2*ZEInv)
            x101= PAQPx1
            x011= QCPQx1
            x111= PAQPx1*QCPQx1 + B001
            x102= PAQPx2
            x012= QCPQx2
            x112= PAQPx2*QCPQx2 + B002
            y101= PAQPy1
            y011= QCPQy1
            y111= PAQPy1*QCPQy1 + B001
            y102= PAQPy2
            y012= QCPQy2
            y112= PAQPy2*QCPQy2 + B002
            z101= PAQPz1
            z011= QCPQz1
            z111= PAQPz1*QCPQz1 + B001
            z102= PAQPz2
            z012= QCPQz2
            z112= PAQPz2*QCPQz2 + B002
            PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv)
            EFInt(iZeta,iEta,1,1) = PreFct * (x111     *w1+x112     *w2)
            EFInt(iZeta,iEta,1,2) = PreFct * (x101*y011*w1+x102*y012*w2)
            EFInt(iZeta,iEta,1,3) = PreFct * (x101*z011*w1+x102*z012*w2)
            EFInt(iZeta,iEta,2,1) = PreFct * (y101*x011*w1+y102*x012*w2)
            EFInt(iZeta,iEta,2,2) = PreFct * (y111     *w1 + y112   *w2)
            EFInt(iZeta,iEta,2,3) = PreFct * (y101*z011*w1+y102*z012*w2)
            EFInt(iZeta,iEta,3,1) = PreFct * (z101*x011*w1+z102*x012*w2)
            EFInt(iZeta,iEta,3,2) = PreFct * (z101*y011*w1+z102*y012*w2)
            EFInt(iZeta,iEta,3,3) = PreFct * (z111     *w1 + z112   *w2)
 20      Continue
 10   Continue
      Go To 99
*
*-----AAAA case
*
 100  Continue
      z = - x0(1)
      w1=(((((CW6(1,1)*z+CW5(1,1))*z+CW4(1,1))*z+CW3(1,1))*z+
     &   CW2(1,1))*z+CW1(1,1))*z+Cw0(1,1)
      w2=(((((CW6(1,2)*z+CW5(1,2))*z+CW4(1,2))*z+CW3(1,2))*z+
     &   CW2(1,2))*z+CW1(1,2))*z+Cw0(1,2)
      r1=(((((CR6(1,1)*z+CR5(1,1))*z+CR4(1,1))*z+CR3(1,1))*z+
     &   CR2(1,1))*z+CR1(1,1))*z+CR0(1,1)
      r2=(((((CR6(1,2)*z+CR5(1,2))*z+CR4(1,2))*z+CR3(1,2))*z+
     &   CR2(1,2))*z+CR1(1,2))*z+CR0(1,2)
      Do 11 iEta = 1, nEta
         Do 21 iZeta = 1, nZeta
            ZEInv = One/(Eta(iEta)+Zeta(iZeta)
     >              +(Eta(iEta)*Zeta(iZeta)*ChiI2)*Dble(IsChi))
            B001 = w1* Half * (r1*ZEInv)
            B002 = w2* Half * (r2*ZEInv)
            PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv)
            EFInt(iZeta,iEta,1,1) = (PreFct * (B001+B002))
            EFInt(iZeta,iEta,1,2) = Zero
            EFInt(iZeta,iEta,1,3) = Zero
            EFInt(iZeta,iEta,2,1) = Zero
            EFInt(iZeta,iEta,2,2) = (PreFct * (B001+B002))
            EFInt(iZeta,iEta,2,3) = Zero
            EFInt(iZeta,iEta,3,1) = Zero
            EFInt(iZeta,iEta,3,2) = Zero
            EFInt(iZeta,iEta,3,3) = (PreFct * (B001+B002))
 21      Continue
 11   Continue
      Go To 99
*
*-----AACD case
*
 200  Continue
      Do 12 iEta = 1, nEta
         Do 22 iZeta = 1, nZeta
         PQx = CoorAC(1,1)-Q(iEta,1)
         PQy = CoorAC(2,1)-Q(iEta,2)
         PQz = CoorAC(3,1)-Q(iEta,3)
            ZEInv = One/(Eta(iEta)+Zeta(iZeta)
     >              +(Eta(iEta)*Zeta(iZeta)*ChiI2)*Dble(IsChi))
            rho = Zeta(iZeta)*Eta(iEta)*ZEInv
            T = rho * (PQx**2 + PQy**2 + PQz**2)
            If (T.lt.TMax) Then
               n = iPntr(Int((T+dddx)*xdInv))
               z = T - x0(n)
               w1=(((((CW6(n,1)*z+CW5(n,1))*z+CW4(n,1))*z+CW3(n,1))*z+
     &            CW2(n,1))*z+CW1(n,1))*z+Cw0(n,1)
               w2=(((((CW6(n,2)*z+CW5(n,2))*z+CW4(n,2))*z+CW3(n,2))*z+
     &            CW2(n,2))*z+CW1(n,2))*z+Cw0(n,2)
               r1=(((((CR6(n,1)*z+CR5(n,1))*z+CR4(n,1))*z+CR3(n,1))*z+
     &            CR2(n,1))*z+CR1(n,1))*z+CR0(n,1)
               r2=(((((CR6(n,2)*z+CR5(n,2))*z+CR4(n,2))*z+CR3(n,2))*z+
     &            CR2(n,2))*z+CR1(n,2))*z+CR0(n,2)
            Else
               ai = 1.0D0/T
               si = Sqrt(ai)
               w1= HerW(1)*si
               w2= HerW(2)*si
               r1= HerR2(1)*ai
               r2= HerR2(2)*ai
            End If
            Eu21 =   Eta(iEta) * (r1*ZEInv)
            Eu22 =   Eta(iEta) * (r2*ZEInv)
            Zu21 = Zeta(iZeta) * (r1*ZEInv)
            Zu22 = Zeta(iZeta) * (r2*ZEInv)
            PAQPx1 = - Eu21 * PQx
            PAQPx2 = - Eu22 * PQx
            PAQPy1 = - Eu21 * PQy
            PAQPy2 = - Eu22 * PQy
            PAQPz1 = - Eu21 * PQz
            PAQPz2 = - Eu22 * PQz
            QCPQx1 = (Q(iEta,1)  - CoorAC(1,2)) + Zu21 * PQx
            QCPQx2 = (Q(iEta,1)  - CoorAC(1,2)) + Zu22 * PQx
            QCPQy1 = (Q(iEta,2)  - CoorAC(2,2)) + Zu21 * PQy
            QCPQy2 = (Q(iEta,2)  - CoorAC(2,2)) + Zu22 * PQy
            QCPQz1 = (Q(iEta,3)  - CoorAC(3,2)) + Zu21 * PQz
            QCPQz2 = (Q(iEta,3)  - CoorAC(3,2)) + Zu22 * PQz
            B001 = Half * (r1*ZEInv)
            B002 = Half * (r2*ZEInv)
            x101= PAQPx1
            x011= QCPQx1
            x111= PAQPx1*QCPQx1 + B001
            x102= PAQPx2
            x012= QCPQx2
            x112= PAQPx2*QCPQx2 + B002
            y101= PAQPy1
            y011= QCPQy1
            y111= PAQPy1*QCPQy1 + B001
            y102= PAQPy2
            y012= QCPQy2
            y112= PAQPy2*QCPQy2 + B002
            z101= PAQPz1
            z011= QCPQz1
            z111= PAQPz1*QCPQz1 + B001
            z102= PAQPz2
            z012= QCPQz2
            z112= PAQPz2*QCPQz2 + B002
            PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv)
            EFInt(iZeta,iEta,1,1) = PreFct * (x111     *w1+x112     *w2)
            EFInt(iZeta,iEta,1,2) = PreFct * (x101*y011*w1+x102*y012*w2)
            EFInt(iZeta,iEta,1,3) = PreFct * (x101*z011*w1+x102*z012*w2)
            EFInt(iZeta,iEta,2,1) = PreFct * (y101*x011*w1+y102*x012*w2)
            EFInt(iZeta,iEta,2,2) = PreFct * (y111     *w1 + y112   *w2)
            EFInt(iZeta,iEta,2,3) = PreFct * (y101*z011*w1+y102*z012*w2)
            EFInt(iZeta,iEta,3,1) = PreFct * (z101*x011*w1+z102*x012*w2)
            EFInt(iZeta,iEta,3,2) = PreFct * (z101*y011*w1+z102*y012*w2)
            EFInt(iZeta,iEta,3,3) = PreFct * (z111     *w1 + z112   *w2)
 22      Continue
 12   Continue
      Go To 99
*
*-----ABCC case
*
 300  Continue
      Do 13 iEta = 1, nEta
         Do 23 iZeta = 1, nZeta
            ZEInv = One/(Eta(iEta)+Zeta(iZeta)
     >              +(Eta(iEta)*Zeta(iZeta)*ChiI2)*Dble(IsChi))
            rho = Zeta(iZeta)*Eta(iEta)*ZEInv
            PQx = P(iZeta,1)-CoorAC(1,2)
            PQy = P(iZeta,2)-CoorAC(2,2)
            PQz = P(iZeta,3)-CoorAC(3,2)
            T = rho * (PQx**2 + PQy**2 + PQz**2)
            If (T.lt.TMax) Then
               n = iPntr(Int((T+dddx)*xdInv))
               z = T - x0(n)
               w1=(((((CW6(n,1)*z+CW5(n,1))*z+CW4(n,1))*z+CW3(n,1))*z+
     &            CW2(n,1))*z+CW1(n,1))*z+Cw0(n,1)
               w2=(((((CW6(n,2)*z+CW5(n,2))*z+CW4(n,2))*z+CW3(n,2))*z+
     &            CW2(n,2))*z+CW1(n,2))*z+Cw0(n,2)
               r1=(((((CR6(n,1)*z+CR5(n,1))*z+CR4(n,1))*z+CR3(n,1))*z+
     &            CR2(n,1))*z+CR1(n,1))*z+CR0(n,1)
               r2=(((((CR6(n,2)*z+CR5(n,2))*z+CR4(n,2))*z+CR3(n,2))*z+
     &            CR2(n,2))*z+CR1(n,2))*z+CR0(n,2)
            Else
               ai = 1.0D0/T
               si = Sqrt(ai)
               w1= HerW(1)*si
               w2= HerW(2)*si
               r1= HerR2(1)*ai
               r2= HerR2(2)*ai
            End If
            Eu21 =   Eta(iEta) * (r1*ZEInv)
            Eu22 =   Eta(iEta) * (r2*ZEInv)
            Zu21 = Zeta(iZeta) * (r1*ZEInv)
            Zu22 = Zeta(iZeta) * (r2*ZEInv)
            PAQPx1 = (P(iZeta,1) - CoorAC(1,1)) - Eu21 * PQx
            PAQPx2 = (P(iZeta,1) - CoorAC(1,1)) - Eu22 * PQx
            PAQPy1 = (P(iZeta,2) - CoorAC(2,1)) - Eu21 * PQy
            PAQPy2 = (P(iZeta,2) - CoorAC(2,1)) - Eu22 * PQy
            PAQPz1 = (P(iZeta,3) - CoorAC(3,1)) - Eu21 * PQz
            PAQPz2 = (P(iZeta,3) - CoorAC(3,1)) - Eu22 * PQz
            QCPQx1 = Zu21 * PQx
            QCPQx2 = Zu22 * PQx
            QCPQy1 = Zu21 * PQy
            QCPQy2 = Zu22 * PQy
            QCPQz1 = Zu21 * PQz
            QCPQz2 = Zu22 * PQz
            B001 = Half * (r1*ZEInv)
            B002 = Half * (r2*ZEInv)
            x101= PAQPx1
            x011= QCPQx1
            x111= PAQPx1*QCPQx1 + B001
            x102= PAQPx2
            x012= QCPQx2
            x112= PAQPx2*QCPQx2 + B002
            y101= PAQPy1
            y011= QCPQy1
            y111= PAQPy1*QCPQy1 + B001
            y102= PAQPy2
            y012= QCPQy2
            y112= PAQPy2*QCPQy2 + B002
            z101= PAQPz1
            z011= QCPQz1
            z111= PAQPz1*QCPQz1 + B001
            z102= PAQPz2
            z012= QCPQz2
            z112= PAQPz2*QCPQz2 + B002
            PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv)
            EFInt(iZeta,iEta,1,1) = PreFct * (x111     *w1+x112     *w2)
            EFInt(iZeta,iEta,1,2) = PreFct * (x101*y011*w1+x102*y012*w2)
            EFInt(iZeta,iEta,1,3) = PreFct * (x101*z011*w1+x102*z012*w2)
            EFInt(iZeta,iEta,2,1) = PreFct * (y101*x011*w1+y102*x012*w2)
            EFInt(iZeta,iEta,2,2) = PreFct * (y111     *w1 + y112   *w2)
            EFInt(iZeta,iEta,2,3) = PreFct * (y101*z011*w1+y102*z012*w2)
            EFInt(iZeta,iEta,3,1) = PreFct * (z101*x011*w1+z102*x012*w2)
            EFInt(iZeta,iEta,3,2) = PreFct * (z101*y011*w1+z102*y012*w2)
            EFInt(iZeta,iEta,3,3) = PreFct * (z111     *w1 + z112   *w2)
 23      Continue
 13   Continue
      Go To 99
*
*-----AACC case
*
 400  Continue
      PQx = CoorAC(1,1)-CoorAC(1,2)
      PQy = CoorAC(2,1)-CoorAC(2,2)
      PQz = CoorAC(3,1)-CoorAC(3,2)
      PQ2 = (PQx**2 + PQy**2 + PQz**2)
      Do 14 iEta = 1, nEta
         Do 24 iZeta = 1, nZeta
            ZEInv = One/(Eta(iEta)+Zeta(iZeta)
     >              +(Eta(iEta)*Zeta(iZeta)*ChiI2)*Dble(IsChi))
            rho = Zeta(iZeta)*Eta(iEta)*ZEInv
            T = rho * PQ2
            If (T.lt.TMax) Then
               n = iPntr(Int((T+dddx)*xdInv))
               z = T - x0(n)
               w1=(((((CW6(n,1)*z+CW5(n,1))*z+CW4(n,1))*z+CW3(n,1))*z+
     &            CW2(n,1))*z+CW1(n,1))*z+Cw0(n,1)
               w2=(((((CW6(n,2)*z+CW5(n,2))*z+CW4(n,2))*z+CW3(n,2))*z+
     &            CW2(n,2))*z+CW1(n,2))*z+Cw0(n,2)
               r1=(((((CR6(n,1)*z+CR5(n,1))*z+CR4(n,1))*z+CR3(n,1))*z+
     &            CR2(n,1))*z+CR1(n,1))*z+CR0(n,1)
               r2=(((((CR6(n,2)*z+CR5(n,2))*z+CR4(n,2))*z+CR3(n,2))*z+
     &            CR2(n,2))*z+CR1(n,2))*z+CR0(n,2)
            Else
               ai = 1.0D0/T
               si = Sqrt(ai)
               w1= HerW(1)*si
               w2= HerW(2)*si
               r1= HerR2(1)*ai
               r2= HerR2(2)*ai
            End If
            Eu21 =   Eta(iEta) * (r1*ZEInv)
            Eu22 =   Eta(iEta) * (r2*ZEInv)
            Zu21 = Zeta(iZeta) * (r1*ZEInv)
            Zu22 = Zeta(iZeta) * (r2*ZEInv)
            PAQPx1 = - Eu21 * PQx
            PAQPx2 = - Eu22 * PQx
            PAQPy1 = - Eu21 * PQy
            PAQPy2 = - Eu22 * PQy
            PAQPz1 = - Eu21 * PQz
            PAQPz2 = - Eu22 * PQz
            QCPQx1 =   Zu21 * PQx
            QCPQx2 =   Zu22 * PQx
            QCPQy1 =   Zu21 * PQy
            QCPQy2 =   Zu22 * PQy
            QCPQz1 =   Zu21 * PQz
            QCPQz2 =   Zu22 * PQz
            B001 = Half * (r1*ZEInv)
            B002 = Half * (r2*ZEInv)
            x101= PAQPx1
            x011= QCPQx1
            x111= PAQPx1*QCPQx1 + B001
            x102= PAQPx2
            x012= QCPQx2
            x112= PAQPx2*QCPQx2 + B002
            y101= PAQPy1
            y011= QCPQy1
            y111= PAQPy1*QCPQy1 + B001
            y102= PAQPy2
            y012= QCPQy2
            y112= PAQPy2*QCPQy2 + B002
            z101= PAQPz1
            z011= QCPQz1
            z111= PAQPz1*QCPQz1 + B001
            z102= PAQPz2
            z012= QCPQz2
            z112= PAQPz2*QCPQz2 + B002
            PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv)
            EFInt(iZeta,iEta,1,1) = PreFct * (x111     *w1+x112     *w2)
            EFInt(iZeta,iEta,1,2) = PreFct * (x101*y011*w1+x102*y012*w2)
            EFInt(iZeta,iEta,1,3) = PreFct * (x101*z011*w1+x102*z012*w2)
            EFInt(iZeta,iEta,2,1) = PreFct * (y101*x011*w1+y102*x012*w2)
            EFInt(iZeta,iEta,2,2) = PreFct * (y111     *w1 + y112   *w2)
            EFInt(iZeta,iEta,2,3) = PreFct * (y101*z011*w1+y102*z012*w2)
            EFInt(iZeta,iEta,3,1) = PreFct * (z101*x011*w1+z102*x012*w2)
            EFInt(iZeta,iEta,3,2) = PreFct * (z101*y011*w1+z102*y012*w2)
            EFInt(iZeta,iEta,3,3) = PreFct * (z111     *w1 + z112   *w2)
 24      Continue
 14   Continue
      Go To 99
*
 99   Continue
*
*     Call qExit('psps')
      Return
      End
