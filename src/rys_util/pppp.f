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
      Subroutine pppp(EFInt,Zeta,ZInv,nZeta,P,lP,rKappAB,A,B,
     &                       Eta,EInv, nEta,Q,lQ,rKappCD,C,D,
     &                CoorAC,TMax,
     &                iPntr,nPntr,x0,nMax,CW6,CW5,CW4,CW3,CW2,CW1,CW0,
     &                                    CR6,CR5,CR4,CR3,CR2,CR1,CR0,
     &                ddx,HerW,HerR2,IsChi,ChiI2)
************************************************************************
*                                                                      *
* Object: to compute the primitive integrals of type (pp|pp).          *
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
      Real*8 EFInt(nZeta,nEta,81), Zeta(nZeta), Eta(nEta),
     &       CoorAC(3,2), ZInv(nZeta), EInv(nEta),
     &       P(lP,3), Q(lQ,3), A(3), B(3), C(3), D(3),
     &       rKappAB(nZeta), rKappCD(nEta),
     &       x0(nMax),
     &       CW6(nMax,3), CW5(nMax,3), CW4(nMax,3), CW3(nMax,3),
     &       CW2(nMax,3), CW1(nMax,3), CW0(nMax,3),
     &       CR6(nMax,3), CR5(nMax,3), CR4(nMax,3), CR3(nMax,3),
     &       CR2(nMax,3), CR1(nMax,3), CR0(nMax,3),
     &       HerW(3), HerR2(3)
      Integer iPntr(nPntr)
      Logical ABeqCD, EQ
*
*     Call qEnter('pppp')
*
      xdInv=One/ddx
      dddx = ddx/10d0 + ddx
*
      ABeqCD = EQ(A,B) .and. EQ(A,C) .and. EQ(A,D)
      If (            ABeqCD           ) Go To 100
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
            rho = Zeta(iZeta)*(Eta(iEta)*ZEInv)
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
               w3=(((((CW6(n,3)*z+CW5(n,3))*z+CW4(n,3))*z+CW3(n,3))*z+
     &            CW2(n,3))*z+CW1(n,3))*z+Cw0(n,3)
               r1=(((((CR6(n,1)*z+CR5(n,1))*z+CR4(n,1))*z+CR3(n,1))*z+
     &            CR2(n,1))*z+CR1(n,1))*z+CR0(n,1)
               r2=(((((CR6(n,2)*z+CR5(n,2))*z+CR4(n,2))*z+CR3(n,2))*z+
     &            CR2(n,2))*z+CR1(n,2))*z+CR0(n,2)
               r3=(((((CR6(n,3)*z+CR5(n,3))*z+CR4(n,3))*z+CR3(n,3))*z+
     &            CR2(n,3))*z+CR1(n,3))*z+CR0(n,3)
            Else
               ai = 1.0D0/T
               si = Sqrt(ai)
               r1= HerR2(1)*ai
               r2= HerR2(2)*ai
               r3= HerR2(3)*ai
               w1= HerW(1)*si
               w2= HerW(2)*si
               w3= HerW(3)*si
            End If
            PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv)
            w1 = PreFct*w1
            w2 = PreFct*w2
            w3 = PreFct*w3
            Eu21 =     Eta(iEta)*(r1*ZEInv)
            Eu22 =     Eta(iEta)*(r2*ZEInv)
            Eu23 =     Eta(iEta)*(r3*ZEInv)
            Zu21 =   Zeta(iZeta)*(r1*ZEInv)
            Zu22 =   Zeta(iZeta)*(r2*ZEInv)
            Zu23 =   Zeta(iZeta)*(r3*ZEInv)
            PAQPx1 = (P(iZeta,1) - CoorAC(1,1)) - Eu21 * PQx
            PAQPx2 = (P(iZeta,1) - CoorAC(1,1)) - Eu22 * PQx
            PAQPx3 = (P(iZeta,1) - CoorAC(1,1)) - Eu23 * PQx
            PAQPy1 = (P(iZeta,2) - CoorAC(2,1)) - Eu21 * PQy
            PAQPy2 = (P(iZeta,2) - CoorAC(2,1)) - Eu22 * PQy
            PAQPy3 = (P(iZeta,2) - CoorAC(2,1)) - Eu23 * PQy
            PAQPz1 = (P(iZeta,3) - CoorAC(3,1)) - Eu21 * PQz
            PAQPz2 = (P(iZeta,3) - CoorAC(3,1)) - Eu22 * PQz
            PAQPz3 = (P(iZeta,3) - CoorAC(3,1)) - Eu23 * PQz
            QCPQx1 = ( Q(iEta,1) - CoorAC(1,2)) + Zu21 * PQx
            QCPQx2 = ( Q(iEta,1) - CoorAC(1,2)) + Zu22 * PQx
            QCPQx3 = ( Q(iEta,1) - CoorAC(1,2)) + Zu23 * PQx
            QCPQy1 = ( Q(iEta,2) - CoorAC(2,2)) + Zu21 * PQy
            QCPQy2 = ( Q(iEta,2) - CoorAC(2,2)) + Zu22 * PQy
            QCPQy3 = ( Q(iEta,2) - CoorAC(2,2)) + Zu23 * PQy
            QCPQz1 = ( Q(iEta,3) - CoorAC(3,2)) + Zu21 * PQz
            QCPQz2 = ( Q(iEta,3) - CoorAC(3,2)) + Zu22 * PQz
            QCPQz3 = ( Q(iEta,3) - CoorAC(3,2)) + Zu23 * PQz
            B001 = Half * (r1*ZEInv)
            B002 = Half * (r2*ZEInv)
            B003 = Half * (r3*ZEInv)
            B101 = (Half - Half * Eu21) * ZInv(iZeta)
            B102 = (Half - Half * Eu22) * ZInv(iZeta)
            B103 = (Half - Half * Eu23) * ZInv(iZeta)
            B011 = (Half - Half * Zu21) * EInv(iEta)
            B012 = (Half - Half * Zu22) * EInv(iEta)
            B013 = (Half - Half * Zu23) * EInv(iEta)
            x101= PAQPx1
            x102= PAQPx2
            x103= PAQPx3
            x011= QCPQx1
            x012= QCPQx2
            x013= QCPQx3
            x201= PAQPx1*x101 + B101
            x202= PAQPx2*x102 + B102
            x203= PAQPx3*x103 + B103
            x021= QCPQx1*x011 + B011
            x022= QCPQx2*x012 + B012
            x023= QCPQx3*x013 + B013
            x111= PAQPx1*QCPQx1 + B001
            x112= PAQPx2*QCPQx2 + B002
            x113= PAQPx3*QCPQx3 + B003
            x211= PAQPx1*x111 + B101*x011 + B001*x101
            x212= PAQPx2*x112 + B102*x012 + B002*x102
            x213= PAQPx3*x113 + B103*x013 + B003*x103
            x121= QCPQx1*x111 + B011*x101 + B001*x011
            x122= QCPQx2*x112 + B012*x102 + B002*x012
            x123= QCPQx3*x113 + B013*x103 + B003*x013
            x221= PAQPx1*x121 + B101*x021 + Two*B001*x111
            x222= PAQPx2*x122 + B102*x022 + Two*B002*x112
            x223= PAQPx3*x123 + B103*x023 + Two*B003*x113
            y101= PAQPy1
            y102= PAQPy2
            y103= PAQPy3
            y011= QCPQy1
            y012= QCPQy2
            y013= QCPQy3
            y201= PAQPy1*y101 + B101
            y202= PAQPy2*y102 + B102
            y203= PAQPy3*y103 + B103
            y021= QCPQy1*y011 + B011
            y022= QCPQy2*y012 + B012
            y023= QCPQy3*y013 + B013
            y111= PAQPy1*QCPQy1 + B001
            y112= PAQPy2*QCPQy2 + B002
            y113= PAQPy3*QCPQy3 + B003
            y211= PAQPy1*y111 + B101*y011 + B001*y101
            y212= PAQPy2*y112 + B102*y012 + B002*y102
            y213= PAQPy3*y113 + B103*y013 + B003*y103
            y121= QCPQy1*y111 + B011*y101 + B001*y011
            y122= QCPQy2*y112 + B012*y102 + B002*y012
            y123= QCPQy3*y113 + B013*y103 + B003*y013
            y221= PAQPy1*y121 + B101*y021 + Two*B001*y111
            y222= PAQPy2*y122 + B102*y022 + Two*B002*y112
            y223= PAQPy3*y123 + B103*y023 + Two*B003*y113
            z101= PAQPz1*w1
            z102= PAQPz2*w2
            z103= PAQPz3*w3
            z011= QCPQz1*w1
            z012= QCPQz2*w2
            z013= QCPQz3*w3
            z201= PAQPz1*z101 + B101*w1
            z202= PAQPz2*z102 + B102*w2
            z203= PAQPz3*z103 + B103*w3
            z021= QCPQz1*z011 + B011*w1
            z022= QCPQz2*z012 + B012*w2
            z023= QCPQz3*z013 + B013*w3
            z111= PAQPz1*QCPQz1*w1 + B001*w1
            z112= PAQPz2*QCPQz2*w2 + B002*w2
            z113= PAQPz3*QCPQz3*w3 + B003*w3
            z211= PAQPz1*z111 + B101*z011 + B001*z101
            z212= PAQPz2*z112 + B102*z012 + B002*z102
            z213= PAQPz3*z113 + B103*z013 + B003*z103
            z121= QCPQz1*z111 + B011*z101 + B001*z011
            z122= QCPQz2*z112 + B012*z102 + B002*z012
            z123= QCPQz3*z113 + B013*z103 + B003*z013
            z221= PAQPz1*z121 + B101*z021 + Two*B001*z111
            z222= PAQPz2*z122 + B102*z022 + Two*B002*z112
            z223= PAQPz3*z123 + B103*z023 + Two*B003*z113
            EFInt(iZeta,iEta, 1)=
     &    (x111        *   w1)+(x112        *   w2)+(x113        *   w3)
            EFInt(iZeta,iEta, 2)=
     &    (x011 * y101)*   w1 +(x012 * y102)*   w2 +(x013 * y103)*   w3
            EFInt(iZeta,iEta, 3)=
     &     x011        * z101 + x012        * z102 + x013        * z103
            EFInt(iZeta,iEta, 4)=
     &    (x211        *   w1)+(x212        *   w2)+(x213        *   w3)
            EFInt(iZeta,iEta, 5)=
     &    (x111 *   w1)* y101 +(x112 *   w2)* y102 +(x113 *   w3)* y103
            EFInt(iZeta,iEta, 6)=
     &     x111        * z101 + x112        * z102 + x113        * z103
            EFInt(iZeta,iEta, 7)=
     &     x011 * y201 *   w1 + x012 * y202 *   w2 + x013 * y203 *   w3
            EFInt(iZeta,iEta, 8)=
     &    (x011 * y101)* z101 +(x012 * y102)* z102 +(x013 * y103)* z103
            EFInt(iZeta,iEta, 9)=
     &     x011        * z201 + x012        * z202 + x013        * z203
            EFInt(iZeta,iEta,10)=
     &    (x101 * y011)*   w1 +(x102 * y012)*   w2 +(x103 * y013)*   w3
            EFInt(iZeta,iEta,11)=
     &            y111 *   w1 +        y112 *   w2 +        y113 *   w3
            EFInt(iZeta,iEta,12)=
     &           (y011 * z101)+       (y012 * z102)+       (y013 * z103)
            EFInt(iZeta,iEta,13)=
     &    (x201 * y011)*   w1 +(x202 * y012)*   w2 +(x203 * y013)*   w3
            EFInt(iZeta,iEta,14)=
     &    (x101 * y111)*   w1 +(x102 * y112)*   w2 +(x103 * y113)*   w3
            EFInt(iZeta,iEta,15)=
     &    (x101 * y011)* z101 +(x102 * y012)* z102 +(x103 * y013)* z103
            EFInt(iZeta,iEta,16)=
     &           (y211 *   w1)+       (y212 *   w2)+       (y213 *   w3)
            EFInt(iZeta,iEta,17)=
     &            y111 * z101 +        y112 * z102 +        y113 * z103
            EFInt(iZeta,iEta,18)=
     &            y011 * z201 +        y012 * z202 +        y013 * z203
            EFInt(iZeta,iEta,19)=
     &     x101        * z011 + x102        * z012 + x103        * z013
            EFInt(iZeta,iEta,20)=
     &           (y101 * z011)+       (y102 * z012)+       (y103 * z013)
            EFInt(iZeta,iEta,21)=
     &                   z111 +               z112 +               z113
            EFInt(iZeta,iEta,22)=
     &     x201        * z011 + x202        * z012 + x203        * z013
            EFInt(iZeta,iEta,23)=
     &     x101 *(y101 * z011)+ x102 *(y102 * z012)+ x103 *(y103 * z013)
            EFInt(iZeta,iEta,24)=
     &     x101        * z111 + x102        * z112 + x103        * z113
            EFInt(iZeta,iEta,25)=
     &            y201 * z011 +        y202 * z012 +        y203 * z013
            EFInt(iZeta,iEta,26)=
     &            y101 * z111 +        y102 * z112 +        y103 * z113
            EFInt(iZeta,iEta,27)=
     &                   z211 +               z212 +               z213
            EFInt(iZeta,iEta,28)=
     &    (x121        *   w1)+(x122        *   w2)+(x123        *   w3)
            EFInt(iZeta,iEta,29)=
     &    (x021 * y101)*   w1 +(x022 * y102)*   w2 +(x023 * y103)*   w3
            EFInt(iZeta,iEta,30)=
     &     x021        * z101 + x022        * z102 + x023        * z103
            EFInt(iZeta,iEta,31)=
     &     x221        *   w1 + x222        *   w2 + x223        *   w3
            EFInt(iZeta,iEta,32)=
     &    (x121 *   w1)* y101 +(x122 *   w2)* y102 +(x123 *   w3)* y103
            EFInt(iZeta,iEta,33)=
     &     x121        * z101 + x122        * z102 + x123        * z103
            EFInt(iZeta,iEta,34)=
     &     x021 * y201 *   w1 + x022 * y202 *   w2 + x023 * y203 *   w3
            EFInt(iZeta,iEta,35)=
     &    (x021 * y101)* z101 +(x022 * y102)* z102 +(x023 * y103)* z103
            EFInt(iZeta,iEta,36)=
     &     x021        * z201 + x022        * z202 + x023        * z203
            EFInt(iZeta,iEta,37)=
     &    (x111 *   w1)* y011 +(x112 *   w2)* y012 +(x113 *   w3)* y013
            EFInt(iZeta,iEta,38)=
     &    (x011 * y111)*   w1 +(x012 * y112)*   w2 +(x013 * y113)*   w3
            EFInt(iZeta,iEta,39)=
     &    (x011 * y011)* z101 +(x012 * y012)* z102 +(x013 * y013)* z103
            EFInt(iZeta,iEta,40)=
     &    (x211 *   w1)* y011 +(x212 *   w2)* y012 +(x213 *   w3)* y013
            EFInt(iZeta,iEta,41)=
     &     x111 * y111 *   w1 + x112 * y112 *   w2 + x113 * y113 *   w3
            EFInt(iZeta,iEta,42)=
     &     x111 *(y011 * z101)+ x112 *(y012 * z102)+ x113 *(y013 * z103)
            EFInt(iZeta,iEta,43)=
     &     x011 *(y211 *   w1)+ x012 *(y212 *   w2)+ x013 *(y213 *   w3)
            EFInt(iZeta,iEta,44)=
     &    (x011 * y111)* z101 +(x012 * y112)* z102 +(x013 * y113)* z103
            EFInt(iZeta,iEta,45)=
     &    (x011 * y011)* z201 +(x012 * y012)* z202 +(x013 * y013)* z203
            EFInt(iZeta,iEta,46)=
     &     x111        * z011 + x112        * z012 + x113        * z013
            EFInt(iZeta,iEta,47)=
     &    (x011 * y101)* z011 +(x012 * y102)* z012 +(x013 * y103)* z013
            EFInt(iZeta,iEta,48)=
     &     x011        * z111 + x012        * z112 + x013        * z113
            EFInt(iZeta,iEta,49)=
     &     x211        * z011 + x212        * z012 + x213        * z013
            EFInt(iZeta,iEta,50)=
     &     x111 *(y101 * z011)+ x112 *(y102 * z012)+ x113 *(y103 * z013)
            EFInt(iZeta,iEta,51)=
     &     x111        * z111 + x112        * z112 + x113        * z113
            EFInt(iZeta,iEta,52)=
     &    (x011 * y201)* z011 +(x012 * y202)* z012 +(x013 * y203)* z013
            EFInt(iZeta,iEta,53)=
     &    (x011 * y101)* z111 +(x012 * y102)* z112 +(x013 * y103)* z113
            EFInt(iZeta,iEta,54)=
     &     x011        * z211 + x012        * z212 + x013        * z213
            EFInt(iZeta,iEta,55)=
     &     x101 * y021 *   w1 + x102 * y022 *   w2 + x103 * y023 *   w3
            EFInt(iZeta,iEta,56)=
     &           (y121 *   w1)+       (y122 *   w2)+       (y123 *   w3)
            EFInt(iZeta,iEta,57)=
     &           (y021 * z101)+       (y022 * z102)+       (y023 * z103)
            EFInt(iZeta,iEta,58)=
     &     x201 * y021 *   w1 + x202 * y022 *   w2 + x203 * y023 *   w3
            EFInt(iZeta,iEta,59)=
     &     x101 *(y121 *   w1)+ x102 *(y122 *   w2)+ x103 *(y123 *   w3)
            EFInt(iZeta,iEta,60)=
     &     x101 *(y021 * z101)+ x102 *(y022 * z102)+ x103 *(y023 * z103)
            EFInt(iZeta,iEta,61)=
     &            y221 *   w1 +        y222 *   w2 +        y223 *   w3
            EFInt(iZeta,iEta,62)=
     &            y121 * z101 +        y122 * z102 +        y123 * z103
            EFInt(iZeta,iEta,63)=
     &            y021 * z201 +        y022 * z202 +        y023 * z203
            EFInt(iZeta,iEta,64)=
     &    (x101 * y011)* z011 +(x102 * y012)* z012 +(x103 * y013)* z013
            EFInt(iZeta,iEta,65)=
     &            y111 * z011 +        y112 * z012 +        y113 * z013
            EFInt(iZeta,iEta,66)=
     &            y011 * z111 +        y012 * z112 +        y013 * z113
            EFInt(iZeta,iEta,67)=
     &    (x201 * y011)* z011 +(x202 * y012)* z012 +(x203 * y013)* z013
            EFInt(iZeta,iEta,68)=
     &    (x101 * y111)* z011 +(x102 * y112)* z012 +(x103 * y113)* z013
            EFInt(iZeta,iEta,69)=
     &    (x101 * y011)* z111 +(x102 * y012)* z112 +(x103 * y013)* z113
            EFInt(iZeta,iEta,70)=
     &            y211 * z011 +        y212 * z012 +        y213 * z013
            EFInt(iZeta,iEta,71)=
     &            y111 * z111 +        y112 * z112 +        y113 * z113
            EFInt(iZeta,iEta,72)=
     &            y011 * z211 +        y012 * z212 +        y013 * z213
            EFInt(iZeta,iEta,73)=
     &     x101        * z021 + x102        * z022 + x103        * z023
            EFInt(iZeta,iEta,74)=
     &           (y101 * z021)+       (y102 * z022)+       (y103 * z023)
            EFInt(iZeta,iEta,75)=
     &                   z121 +               z122 +               z123
            EFInt(iZeta,iEta,76)=
     &     x201        * z021 + x202        * z022 + x203        * z023
            EFInt(iZeta,iEta,77)=
     &     x101 *(y101 * z021)+ x102 *(y102 * z022)+ x103 *(y103 * z023)
            EFInt(iZeta,iEta,78)=
     &     x101        * z121 + x102        * z122 + x103        * z123
            EFInt(iZeta,iEta,79)=
     &            y201 * z021 +        y202 * z022 +        y203 * z023
            EFInt(iZeta,iEta,80)=
     &            y101 * z121 +        y102 * z122 +        y103 * z123
            EFInt(iZeta,iEta,81)=
     &                   z221 +               z222 +               z223
 20      Continue
 10   Continue
      Go To 99
*
*-----AAAA case
*
 100  Continue
      z =  - x0(1)
      ww1=(((((CW6(1,1)*z+CW5(1,1))*z+CW4(1,1))*z+CW3(1,1))*z+
     &   CW2(1,1))*z+CW1(1,1))*z+Cw0(1,1)
      ww2=(((((CW6(1,2)*z+CW5(1,2))*z+CW4(1,2))*z+CW3(1,2))*z+
     &   CW2(1,2))*z+CW1(1,2))*z+Cw0(1,2)
      ww3=(((((CW6(1,3)*z+CW5(1,3))*z+CW4(1,3))*z+CW3(1,3))*z+
     &   CW2(1,3))*z+CW1(1,3))*z+Cw0(1,3)
      r1=(((((CR6(1,1)*z+CR5(1,1))*z+CR4(1,1))*z+CR3(1,1))*z+
     &   CR2(1,1))*z+CR1(1,1))*z+CR0(1,1)
      r2=(((((CR6(1,2)*z+CR5(1,2))*z+CR4(1,2))*z+CR3(1,2))*z+
     &   CR2(1,2))*z+CR1(1,2))*z+CR0(1,2)
      r3=(((((CR6(1,3)*z+CR5(1,3))*z+CR4(1,3))*z+CR3(1,3))*z+
     &   CR2(1,3))*z+CR1(1,3))*z+CR0(1,3)
      Do 11 iEta = 1, nEta
         Do 21 iZeta = 1, nZeta
            ZEInv = One/(Eta(iEta)+Zeta(iZeta)
     >              +(Eta(iEta)*Zeta(iZeta)*ChiI2)*Dble(IsChi))
            PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv)
            w1 = PreFct*ww1
            w2 = PreFct*ww2
            w3 = PreFct*ww3
            Eu21 =     Eta(iEta)*(r1*ZEInv)
            Eu22 =     Eta(iEta)*(r2*ZEInv)
            Eu23 =     Eta(iEta)*(r3*ZEInv)
            Zu21 =   Zeta(iZeta)*(r1*ZEInv)
            Zu22 =   Zeta(iZeta)*(r2*ZEInv)
            Zu23 =   Zeta(iZeta)*(r3*ZEInv)
            B001 = Half * (r1*ZEInv)
            B002 = Half * (r2*ZEInv)
            B003 = Half * (r3*ZEInv)
            B101 = (Half - Half * Eu21) * ZInv(iZeta)
            B102 = (Half - Half * Eu22) * ZInv(iZeta)
            B103 = (Half - Half * Eu23) * ZInv(iZeta)
            B011 = (Half - Half * Zu21) * EInv(iEta)
            B012 = (Half - Half * Zu22) * EInv(iEta)
            B013 = (Half - Half * Zu23) * EInv(iEta)
            x201= B101
            x202= B102
            x203= B103
            x021= B011
            x022= B012
            x023= B013
            x111= B001
            x112= B002
            x113= B003
            x221= B101*x021 + Two*B001*x111
            x222= B102*x022 + Two*B002*x112
            x223= B103*x023 + Two*B003*x113
            y201= B101
            y202= B102
            y203= B103
            y021= B011
            y022= B012
            y023= B013
            y111= B001
            y112= B002
            y113= B003
            y221= B101*y021 + Two*B001*y111
            y222= B102*y022 + Two*B002*y112
            y223= B103*y023 + Two*B003*y113
            z201= B101*w1
            z202= B102*w2
            z203= B103*w3
            z021= B011*w1
            z022= B012*w2
            z023= B013*w3
            z111= B001*w1
            z112= B002*w2
            z113= B003*w3
            z221= B101*z021 + Two*B001*z111
            z222= B102*z022 + Two*B002*z112
            z223= B103*z023 + Two*B003*z113
            EFInt(iZeta,iEta, 1)=
     &     x221        *   w1 + x222        *   w2 + x223        *   w3
            EFInt(iZeta,iEta, 2)= Zero
            EFInt(iZeta,iEta, 3)= Zero
            EFInt(iZeta,iEta, 4)=
     &     x021 * y201 *   w1 + x022 * y202 *   w2 + x023 * y203 *   w3
            EFInt(iZeta,iEta, 5)= Zero
            EFInt(iZeta,iEta, 6)=
     &     x021        * z201 + x022        * z202 + x023        * z203
            EFInt(iZeta,iEta, 7)= Zero
            EFInt(iZeta,iEta, 8)=
     &     x111 * y111 *   w1 + x112 * y112 *   w2 + x113 * y113 *   w3
            EFInt(iZeta,iEta, 9)= Zero
            EFInt(iZeta,iEta,10)= Zero
            EFInt(iZeta,iEta,11)= Zero
            EFInt(iZeta,iEta,12)= Zero
            EFInt(iZeta,iEta,13)= Zero
            EFInt(iZeta,iEta,14)= Zero
            EFInt(iZeta,iEta,15)=
     &     x111        * z111 + x112        * z112 + x113        * z113
            EFInt(iZeta,iEta,16)= Zero
            EFInt(iZeta,iEta,17)= Zero
            EFInt(iZeta,iEta,18)= Zero
            EFInt(iZeta,iEta,19)=
     &     x201 * y021 *   w1 + x202 * y022 *   w2 + x203 * y023 *   w3
            EFInt(iZeta,iEta,20)= Zero
            EFInt(iZeta,iEta,21)= Zero
            EFInt(iZeta,iEta,22)=
     &            y221 *   w1 +        y222 *   w2 +        y223 *   w3
            EFInt(iZeta,iEta,23)= Zero
            EFInt(iZeta,iEta,24)=
     &            y021 * z201 +        y022 * z202 +        y023 * z203
            EFInt(iZeta,iEta,25)= Zero
            EFInt(iZeta,iEta,26)= Zero
            EFInt(iZeta,iEta,27)= Zero
            EFInt(iZeta,iEta,28)= Zero
            EFInt(iZeta,iEta,29)=
     &            y111 * z111 +        y112 * z112 +        y113 * z113
            EFInt(iZeta,iEta,30)= Zero
            EFInt(iZeta,iEta,31)=
     &     x201        * z021 + x202        * z022 + x203        * z023
            EFInt(iZeta,iEta,32)= Zero
            EFInt(iZeta,iEta,33)= Zero
            EFInt(iZeta,iEta,34)=
     &            y201 * z021 +        y202 * z022 +        y203 * z023
            EFInt(iZeta,iEta,35)= Zero
            EFInt(iZeta,iEta,36)=
     &                   z221 +               z222 +               z223
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
         PQ2 = PQx**2 + PQy**2 + PQz**2
            ZEInv = One/(Eta(iEta)+Zeta(iZeta)
     >              +(Eta(iEta)*Zeta(iZeta)*ChiI2)*Dble(IsChi))
            rho = Zeta(iZeta)*(Eta(iEta)*ZEInv)
            T = rho * PQ2
            If (T.lt.TMax) Then
               n = iPntr(Int((T+dddx)*xdInv))
               z = T - x0(n)
               w1=(((((CW6(n,1)*z+CW5(n,1))*z+CW4(n,1))*z+CW3(n,1))*z+
     &            CW2(n,1))*z+CW1(n,1))*z+Cw0(n,1)
               w2=(((((CW6(n,2)*z+CW5(n,2))*z+CW4(n,2))*z+CW3(n,2))*z+
     &            CW2(n,2))*z+CW1(n,2))*z+Cw0(n,2)
               w3=(((((CW6(n,3)*z+CW5(n,3))*z+CW4(n,3))*z+CW3(n,3))*z+
     &            CW2(n,3))*z+CW1(n,3))*z+Cw0(n,3)
               r1=(((((CR6(n,1)*z+CR5(n,1))*z+CR4(n,1))*z+CR3(n,1))*z+
     &            CR2(n,1))*z+CR1(n,1))*z+CR0(n,1)
               r2=(((((CR6(n,2)*z+CR5(n,2))*z+CR4(n,2))*z+CR3(n,2))*z+
     &            CR2(n,2))*z+CR1(n,2))*z+CR0(n,2)
               r3=(((((CR6(n,3)*z+CR5(n,3))*z+CR4(n,3))*z+CR3(n,3))*z+
     &            CR2(n,3))*z+CR1(n,3))*z+CR0(n,3)
            Else
               ai = 1.0D0/T
               si = Sqrt(ai)
               w1= HerW(1)*si
               w2= HerW(2)*si
               w3= HerW(3)*si
               r1= HerR2(1)*ai
               r2= HerR2(2)*ai
               r3= HerR2(3)*ai
            End If
            PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv)
            w1 = PreFct*w1
            w2 = PreFct*w2
            w3 = PreFct*w3
            Eu21 =     Eta(iEta)*(r1*ZEInv)
            Eu22 =     Eta(iEta)*(r2*ZEInv)
            Eu23 =     Eta(iEta)*(r3*ZEInv)
            Zu21 =   Zeta(iZeta)*(r1*ZEInv)
            Zu22 =   Zeta(iZeta)*(r2*ZEInv)
            Zu23 =   Zeta(iZeta)*(r3*ZEInv)
            PAQPx1 = - Eu21 * PQx
            PAQPx2 = - Eu22 * PQx
            PAQPx3 = - Eu23 * PQx
            PAQPy1 = - Eu21 * PQy
            PAQPy2 = - Eu22 * PQy
            PAQPy3 = - Eu23 * PQy
            PAQPz1 = - Eu21 * PQz
            PAQPz2 = - Eu22 * PQz
            PAQPz3 = - Eu23 * PQz
            QCPQx1 = ( Q(iEta,1) - CoorAC(1,2)) + Zu21 * PQx
            QCPQx2 = ( Q(iEta,1) - CoorAC(1,2)) + Zu22 * PQx
            QCPQx3 = ( Q(iEta,1) - CoorAC(1,2)) + Zu23 * PQx
            QCPQy1 = ( Q(iEta,2) - CoorAC(2,2)) + Zu21 * PQy
            QCPQy2 = ( Q(iEta,2) - CoorAC(2,2)) + Zu22 * PQy
            QCPQy3 = ( Q(iEta,2) - CoorAC(2,2)) + Zu23 * PQy
            QCPQz1 = ( Q(iEta,3) - CoorAC(3,2)) + Zu21 * PQz
            QCPQz2 = ( Q(iEta,3) - CoorAC(3,2)) + Zu22 * PQz
            QCPQz3 = ( Q(iEta,3) - CoorAC(3,2)) + Zu23 * PQz
            B001 = Half * (r1*ZEInv)
            B002 = Half * (r2*ZEInv)
            B003 = Half * (r3*ZEInv)
            B101 = (Half - Half * Eu21) * ZInv(iZeta)
            B102 = (Half - Half * Eu22) * ZInv(iZeta)
            B103 = (Half - Half * Eu23) * ZInv(iZeta)
            B011 = (Half - Half * Zu21) * EInv(iEta)
            B012 = (Half - Half * Zu22) * EInv(iEta)
            B013 = (Half - Half * Zu23) * EInv(iEta)
            x101= PAQPx1
            x102= PAQPx2
            x103= PAQPx3
            x011= QCPQx1
            x012= QCPQx2
            x013= QCPQx3
            x201= PAQPx1*x101 + B101
            x202= PAQPx2*x102 + B102
            x203= PAQPx3*x103 + B103
            x021= QCPQx1*x011 + B011
            x022= QCPQx2*x012 + B012
            x023= QCPQx3*x013 + B013
            x111= PAQPx1*QCPQx1 + B001
            x112= PAQPx2*QCPQx2 + B002
            x113= PAQPx3*QCPQx3 + B003
            x211= PAQPx1*x111 + B101*x011 + B001*x101
            x212= PAQPx2*x112 + B102*x012 + B002*x102
            x213= PAQPx3*x113 + B103*x013 + B003*x103
            x121= QCPQx1*x111 + B011*x101 + B001*x011
            x122= QCPQx2*x112 + B012*x102 + B002*x012
            x123= QCPQx3*x113 + B013*x103 + B003*x013
            x221= PAQPx1*x121 + B101*x021 + Two*B001*x111
            x222= PAQPx2*x122 + B102*x022 + Two*B002*x112
            x223= PAQPx3*x123 + B103*x023 + Two*B003*x113
            y101= PAQPy1
            y102= PAQPy2
            y103= PAQPy3
            y011= QCPQy1
            y012= QCPQy2
            y013= QCPQy3
            y201= PAQPy1*y101 + B101
            y202= PAQPy2*y102 + B102
            y203= PAQPy3*y103 + B103
            y021= QCPQy1*y011 + B011
            y022= QCPQy2*y012 + B012
            y023= QCPQy3*y013 + B013
            y111= PAQPy1*QCPQy1 + B001
            y112= PAQPy2*QCPQy2 + B002
            y113= PAQPy3*QCPQy3 + B003
            y211= PAQPy1*y111 + B101*y011 + B001*y101
            y212= PAQPy2*y112 + B102*y012 + B002*y102
            y213= PAQPy3*y113 + B103*y013 + B003*y103
            y121= QCPQy1*y111 + B011*y101 + B001*y011
            y122= QCPQy2*y112 + B012*y102 + B002*y012
            y123= QCPQy3*y113 + B013*y103 + B003*y013
            y221= PAQPy1*y121 + B101*y021 + Two*B001*y111
            y222= PAQPy2*y122 + B102*y022 + Two*B002*y112
            y223= PAQPy3*y123 + B103*y023 + Two*B003*y113
            z101= PAQPz1*w1
            z102= PAQPz2*w2
            z103= PAQPz3*w3
            z011= QCPQz1*w1
            z012= QCPQz2*w2
            z013= QCPQz3*w3
            z201= PAQPz1*z101 + B101*w1
            z202= PAQPz2*z102 + B102*w2
            z203= PAQPz3*z103 + B103*w3
            z021= QCPQz1*z011 + B011*w1
            z022= QCPQz2*z012 + B012*w2
            z023= QCPQz3*z013 + B013*w3
            z111= PAQPz1*QCPQz1*w1 + B001*w1
            z112= PAQPz2*QCPQz2*w2 + B002*w2
            z113= PAQPz3*QCPQz3*w3 + B003*w3
            z211= PAQPz1*z111 + B101*z011 + B001*z101
            z212= PAQPz2*z112 + B102*z012 + B002*z102
            z213= PAQPz3*z113 + B103*z013 + B003*z103
            z121= QCPQz1*z111 + B011*z101 + B001*z011
            z122= QCPQz2*z112 + B012*z102 + B002*z012
            z123= QCPQz3*z113 + B013*z103 + B003*z013
            z221= PAQPz1*z121 + B101*z021 + Two*B001*z111
            z222= PAQPz2*z122 + B102*z022 + Two*B002*z112
            z223= PAQPz3*z123 + B103*z023 + Two*B003*z113
            EFInt(iZeta,iEta, 1)=
     &    (x211        *   w1)+(x212        *   w2)+(x213        *   w3)
            EFInt(iZeta,iEta, 2)=
     &    (x111 *   w1)* y101 +(x112 *   w2)* y102 +(x113 *   w3)* y103
            EFInt(iZeta,iEta, 3)=
     &     x111        * z101 + x112        * z102 + x113        * z103
            EFInt(iZeta,iEta, 4)=
     &     x011 * y201 *   w1 + x012 * y202 *   w2 + x013 * y203 *   w3
            EFInt(iZeta,iEta, 5)=
     &    (x011 * y101)* z101 +(x012 * y102)* z102 +(x013 * y103)* z103
            EFInt(iZeta,iEta, 6)=
     &     x011        * z201 + x012        * z202 + x013        * z203
            EFInt(iZeta,iEta, 7)=
     &    (x201 * y011)*   w1 +(x202 * y012)*   w2 +(x203 * y013)*   w3
            EFInt(iZeta,iEta, 8)=
     &    (x101 * y111)*   w1 +(x102 * y112)*   w2 +(x103 * y113)*   w3
            EFInt(iZeta,iEta, 9)=
     &    (x101 * y011)* z101 +(x102 * y012)* z102 +(x103 * y013)* z103
            EFInt(iZeta,iEta,10)=
     &           (y211 *   w1)+       (y212 *   w2)+       (y213 *   w3)
            EFInt(iZeta,iEta,11)=
     &            y111 * z101 +        y112 * z102 +        y113 * z103
            EFInt(iZeta,iEta,12)=
     &            y011 * z201 +        y012 * z202 +        y013 * z203
            EFInt(iZeta,iEta,13)=
     &     x201        * z011 + x202        * z012 + x203        * z013
            EFInt(iZeta,iEta,14)=
     &     x101 *(y101 * z011)+ x102 *(y102 * z012)+ x103 *(y103 * z013)
            EFInt(iZeta,iEta,15)=
     &     x101        * z111 + x102        * z112 + x103        * z113
            EFInt(iZeta,iEta,16)=
     &            y201 * z011 +        y202 * z012 +        y203 * z013
            EFInt(iZeta,iEta,17)=
     &            y101 * z111 +        y102 * z112 +        y103 * z113
            EFInt(iZeta,iEta,18)=
     &                   z211 +               z212 +               z213
            EFInt(iZeta,iEta,19)=
     &     x221        *   w1 + x222        *   w2 + x223        *   w3
            EFInt(iZeta,iEta,20)=
     &    (x121 *   w1)* y101 +(x122 *   w2)* y102 +(x123 *   w3)* y103
            EFInt(iZeta,iEta,21)=
     &     x121        * z101 + x122        * z102 + x123        * z103
            EFInt(iZeta,iEta,22)=
     &     x021 * y201 *   w1 + x022 * y202 *   w2 + x023 * y203 *   w3
            EFInt(iZeta,iEta,23)=
     &    (x021 * y101)* z101 +(x022 * y102)* z102 +(x023 * y103)* z103
            EFInt(iZeta,iEta,24)=
     &     x021        * z201 + x022        * z202 + x023        * z203
            EFInt(iZeta,iEta,25)=
     &    (x211 *   w1)* y011 +(x212 *   w2)* y012 +(x213 *   w3)* y013
            EFInt(iZeta,iEta,26)=
     &     x111 * y111 *   w1 + x112 * y112 *   w2 + x113 * y113 *   w3
            EFInt(iZeta,iEta,27)=
     &     x111 *(y011 * z101)+ x112 *(y012 * z102)+ x113 *(y013 * z103)
            EFInt(iZeta,iEta,28)=
     &     x011 *(y211 *   w1)+ x012 *(y212 *   w2)+ x013 *(y213 *   w3)
            EFInt(iZeta,iEta,29)=
     &    (x011 * y111)* z101 +(x012 * y112)* z102 +(x013 * y113)* z103
            EFInt(iZeta,iEta,30)=
     &    (x011 * y011)* z201 +(x012 * y012)* z202 +(x013 * y013)* z203
            EFInt(iZeta,iEta,31)=
     &     x211        * z011 + x212        * z012 + x213        * z013
            EFInt(iZeta,iEta,32)=
     &     x111 *(y101 * z011)+ x112 *(y102 * z012)+ x113 *(y103 * z013)
            EFInt(iZeta,iEta,33)=
     &     x111        * z111 + x112        * z112 + x113        * z113
            EFInt(iZeta,iEta,34)=
     &    (x011 * y201)* z011 +(x012 * y202)* z012 +(x013 * y203)* z013
            EFInt(iZeta,iEta,35)=
     &    (x011 * y101)* z111 +(x012 * y102)* z112 +(x013 * y103)* z113
            EFInt(iZeta,iEta,36)=
     &     x011        * z211 + x012        * z212 + x013        * z213
            EFInt(iZeta,iEta,37)=
     &     x201 * y021 *   w1 + x202 * y022 *   w2 + x203 * y023 *   w3
            EFInt(iZeta,iEta,38)=
     &     x101 *(y121 *   w1)+ x102 *(y122 *   w2)+ x103 *(y123 *   w3)
            EFInt(iZeta,iEta,39)=
     &     x101 *(y021 * z101)+ x102 *(y022 * z102)+ x103 *(y023 * z103)
            EFInt(iZeta,iEta,40)=
     &            y221 *   w1 +        y222 *   w2 +        y223 *   w3
            EFInt(iZeta,iEta,41)=
     &            y121 * z101 +        y122 * z102 +        y123 * z103
            EFInt(iZeta,iEta,42)=
     &            y021 * z201 +        y022 * z202 +        y023 * z203
            EFInt(iZeta,iEta,43)=
     &    (x201 * y011)* z011 +(x202 * y012)* z012 +(x203 * y013)* z013
            EFInt(iZeta,iEta,44)=
     &    (x101 * y111)* z011 +(x102 * y112)* z012 +(x103 * y113)* z013
            EFInt(iZeta,iEta,45)=
     &    (x101 * y011)* z111 +(x102 * y012)* z112 +(x103 * y013)* z113
            EFInt(iZeta,iEta,46)=
     &            y211 * z011 +        y212 * z012 +        y213 * z013
            EFInt(iZeta,iEta,47)=
     &            y111 * z111 +        y112 * z112 +        y113 * z113
            EFInt(iZeta,iEta,48)=
     &            y011 * z211 +        y012 * z212 +        y013 * z213
            EFInt(iZeta,iEta,49)=
     &     x201        * z021 + x202        * z022 + x203        * z023
            EFInt(iZeta,iEta,50)=
     &     x101 *(y101 * z021)+ x102 *(y102 * z022)+ x103 *(y103 * z023)
            EFInt(iZeta,iEta,51)=
     &     x101        * z121 + x102        * z122 + x103        * z123
            EFInt(iZeta,iEta,52)=
     &            y201 * z021 +        y202 * z022 +        y203 * z023
            EFInt(iZeta,iEta,53)=
     &            y101 * z121 +        y102 * z122 +        y103 * z123
            EFInt(iZeta,iEta,54)=
     &                   z221 +               z222 +               z223
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
            rho = Zeta(iZeta)*(Eta(iEta)*ZEInv)
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
               w3=(((((CW6(n,3)*z+CW5(n,3))*z+CW4(n,3))*z+CW3(n,3))*z+
     &            CW2(n,3))*z+CW1(n,3))*z+Cw0(n,3)
               r1=(((((CR6(n,1)*z+CR5(n,1))*z+CR4(n,1))*z+CR3(n,1))*z+
     &            CR2(n,1))*z+CR1(n,1))*z+CR0(n,1)
               r2=(((((CR6(n,2)*z+CR5(n,2))*z+CR4(n,2))*z+CR3(n,2))*z+
     &            CR2(n,2))*z+CR1(n,2))*z+CR0(n,2)
               r3=(((((CR6(n,3)*z+CR5(n,3))*z+CR4(n,3))*z+CR3(n,3))*z+
     &            CR2(n,3))*z+CR1(n,3))*z+CR0(n,3)
            Else
               ai = 1.0D0/T
               si = Sqrt(ai)
               w1= HerW(1)*si
               w2= HerW(2)*si
               w3= HerW(3)*si
               r1= HerR2(1)*ai
               r2= HerR2(2)*ai
               r3= HerR2(3)*ai
            End If
            PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv)
            w1 = PreFct*w1
            w2 = PreFct*w2
            w3 = PreFct*w3
            Eu21 =     Eta(iEta)*(r1*ZEInv)
            Eu22 =     Eta(iEta)*(r2*ZEInv)
            Eu23 =     Eta(iEta)*(r3*ZEInv)
            Zu21 =   Zeta(iZeta)*(r1*ZEInv)
            Zu22 =   Zeta(iZeta)*(r2*ZEInv)
            Zu23 =   Zeta(iZeta)*(r3*ZEInv)
            PAQPx1 = (P(iZeta,1) - CoorAC(1,1)) - Eu21 * PQx
            PAQPx2 = (P(iZeta,1) - CoorAC(1,1)) - Eu22 * PQx
            PAQPx3 = (P(iZeta,1) - CoorAC(1,1)) - Eu23 * PQx
            PAQPy1 = (P(iZeta,2) - CoorAC(2,1)) - Eu21 * PQy
            PAQPy2 = (P(iZeta,2) - CoorAC(2,1)) - Eu22 * PQy
            PAQPy3 = (P(iZeta,2) - CoorAC(2,1)) - Eu23 * PQy
            PAQPz1 = (P(iZeta,3) - CoorAC(3,1)) - Eu21 * PQz
            PAQPz2 = (P(iZeta,3) - CoorAC(3,1)) - Eu22 * PQz
            PAQPz3 = (P(iZeta,3) - CoorAC(3,1)) - Eu23 * PQz
            QCPQx1 =  Zu21 * PQx
            QCPQx2 =  Zu22 * PQx
            QCPQx3 =  Zu23 * PQx
            QCPQy1 =  Zu21 * PQy
            QCPQy2 =  Zu22 * PQy
            QCPQy3 =  Zu23 * PQy
            QCPQz1 =  Zu21 * PQz
            QCPQz2 =  Zu22 * PQz
            QCPQz3 =  Zu23 * PQz
            B001 = Half * (r1*ZEInv)
            B002 = Half * (r2*ZEInv)
            B003 = Half * (r3*ZEInv)
            B101 = (Half - Half * Eu21) * ZInv(iZeta)
            B102 = (Half - Half * Eu22) * ZInv(iZeta)
            B103 = (Half - Half * Eu23) * ZInv(iZeta)
            B011 = (Half - Half * Zu21) * EInv(iEta)
            B012 = (Half - Half * Zu22) * EInv(iEta)
            B013 = (Half - Half * Zu23) * EInv(iEta)
            x101= PAQPx1
            x102= PAQPx2
            x103= PAQPx3
            x011= QCPQx1
            x012= QCPQx2
            x013= QCPQx3
            x201= PAQPx1*x101 + B101
            x202= PAQPx2*x102 + B102
            x203= PAQPx3*x103 + B103
            x021= QCPQx1*x011 + B011
            x022= QCPQx2*x012 + B012
            x023= QCPQx3*x013 + B013
            x111= PAQPx1*QCPQx1 + B001
            x112= PAQPx2*QCPQx2 + B002
            x113= PAQPx3*QCPQx3 + B003
            x211= PAQPx1*x111 + B101*x011 + B001*x101
            x212= PAQPx2*x112 + B102*x012 + B002*x102
            x213= PAQPx3*x113 + B103*x013 + B003*x103
            x121= QCPQx1*x111 + B011*x101 + B001*x011
            x122= QCPQx2*x112 + B012*x102 + B002*x012
            x123= QCPQx3*x113 + B013*x103 + B003*x013
            x221= PAQPx1*x121 + B101*x021 + Two*B001*x111
            x222= PAQPx2*x122 + B102*x022 + Two*B002*x112
            x223= PAQPx3*x123 + B103*x023 + Two*B003*x113
            y101= PAQPy1
            y102= PAQPy2
            y103= PAQPy3
            y011= QCPQy1
            y012= QCPQy2
            y013= QCPQy3
            y201= PAQPy1*y101 + B101
            y202= PAQPy2*y102 + B102
            y203= PAQPy3*y103 + B103
            y021= QCPQy1*y011 + B011
            y022= QCPQy2*y012 + B012
            y023= QCPQy3*y013 + B013
            y111= PAQPy1*QCPQy1 + B001
            y112= PAQPy2*QCPQy2 + B002
            y113= PAQPy3*QCPQy3 + B003
            y211= PAQPy1*y111 + B101*y011 + B001*y101
            y212= PAQPy2*y112 + B102*y012 + B002*y102
            y213= PAQPy3*y113 + B103*y013 + B003*y103
            y121= QCPQy1*y111 + B011*y101 + B001*y011
            y122= QCPQy2*y112 + B012*y102 + B002*y012
            y123= QCPQy3*y113 + B013*y103 + B003*y013
            y221= PAQPy1*y121 + B101*y021 + Two*B001*y111
            y222= PAQPy2*y122 + B102*y022 + Two*B002*y112
            y223= PAQPy3*y123 + B103*y023 + Two*B003*y113
            z101= PAQPz1*w1
            z102= PAQPz2*w2
            z103= PAQPz3*w3
            z011= QCPQz1*w1
            z012= QCPQz2*w2
            z013= QCPQz3*w3
            z201= PAQPz1*z101 + B101*w1
            z202= PAQPz2*z102 + B102*w2
            z203= PAQPz3*z103 + B103*w3
            z021= QCPQz1*z011 + B011*w1
            z022= QCPQz2*z012 + B012*w2
            z023= QCPQz3*z013 + B013*w3
            z111= PAQPz1*QCPQz1*w1 + B001*w1
            z112= PAQPz2*QCPQz2*w2 + B002*w2
            z113= PAQPz3*QCPQz3*w3 + B003*w3
            z211= PAQPz1*z111 + B101*z011 + B001*z101
            z212= PAQPz2*z112 + B102*z012 + B002*z102
            z213= PAQPz3*z113 + B103*z013 + B003*z103
            z121= QCPQz1*z111 + B011*z101 + B001*z011
            z122= QCPQz2*z112 + B012*z102 + B002*z012
            z123= QCPQz3*z113 + B013*z103 + B003*z013
            z221= PAQPz1*z121 + B101*z021 + Two*B001*z111
            z222= PAQPz2*z122 + B102*z022 + Two*B002*z112
            z223= PAQPz3*z123 + B103*z023 + Two*B003*z113
            EFInt(iZeta,iEta, 1)=
     &    (x121        *   w1)+(x122        *   w2)+(x123        *   w3)
            EFInt(iZeta,iEta, 2)=
     &    (x021 * y101)*   w1 +(x022 * y102)*   w2 +(x023 * y103)*   w3
            EFInt(iZeta,iEta, 3)=
     &     x021        * z101 + x022        * z102 + x023        * z103
            EFInt(iZeta,iEta, 4)=
     &     x221        *   w1 + x222        *   w2 + x223        *   w3
            EFInt(iZeta,iEta, 5)=
     &    (x121 *   w1)* y101 +(x122 *   w2)* y102 +(x123 *   w3)* y103
            EFInt(iZeta,iEta, 6)=
     &     x121        * z101 + x122        * z102 + x123        * z103
            EFInt(iZeta,iEta, 7)=
     &     x021 * y201 *   w1 + x022 * y202 *   w2 + x023 * y203 *   w3
            EFInt(iZeta,iEta, 8)=
     &    (x021 * y101)* z101 +(x022 * y102)* z102 +(x023 * y103)* z103
            EFInt(iZeta,iEta, 9)=
     &     x021        * z201 + x022        * z202 + x023        * z203
            EFInt(iZeta,iEta,10)=
     &    (x111 *   w1)* y011 +(x112 *   w2)* y012 +(x113 *   w3)* y013
            EFInt(iZeta,iEta,11)=
     &    (x011 * y111)*   w1 +(x012 * y112)*   w2 +(x013 * y113)*   w3
            EFInt(iZeta,iEta,12)=
     &    (x011 * y011)* z101 +(x012 * y012)* z102 +(x013 * y013)* z103
            EFInt(iZeta,iEta,13)=
     &    (x211 *   w1)* y011 +(x212 *   w2)* y012 +(x213 *   w3)* y013
            EFInt(iZeta,iEta,14)=
     &     x111 * y111 *   w1 + x112 * y112 *   w2 + x113 * y113 *   w3
            EFInt(iZeta,iEta,15)=
     &     x111 *(y011 * z101)+ x112 *(y012 * z102)+ x113 *(y013 * z103)
            EFInt(iZeta,iEta,16)=
     &     x011 *(y211 *   w1)+ x012 *(y212 *   w2)+ x013 *(y213 *   w3)
            EFInt(iZeta,iEta,17)=
     &    (x011 * y111)* z101 +(x012 * y112)* z102 +(x013 * y113)* z103
            EFInt(iZeta,iEta,18)=
     &    (x011 * y011)* z201 +(x012 * y012)* z202 +(x013 * y013)* z203
            EFInt(iZeta,iEta,19)=
     &     x111        * z011 + x112        * z012 + x113        * z013
            EFInt(iZeta,iEta,20)=
     &    (x011 * y101)* z011 +(x012 * y102)* z012 +(x013 * y103)* z013
            EFInt(iZeta,iEta,21)=
     &     x011        * z111 + x012        * z112 + x013        * z113
            EFInt(iZeta,iEta,22)=
     &     x211        * z011 + x212        * z012 + x213        * z013
            EFInt(iZeta,iEta,23)=
     &     x111 *(y101 * z011)+ x112 *(y102 * z012)+ x113 *(y103 * z013)
            EFInt(iZeta,iEta,24)=
     &     x111        * z111 + x112        * z112 + x113        * z113
            EFInt(iZeta,iEta,25)=
     &    (x011 * y201)* z011 +(x012 * y202)* z012 +(x013 * y203)* z013
            EFInt(iZeta,iEta,26)=
     &    (x011 * y101)* z111 +(x012 * y102)* z112 +(x013 * y103)* z113
            EFInt(iZeta,iEta,27)=
     &     x011        * z211 + x012        * z212 + x013        * z213
            EFInt(iZeta,iEta,28)=
     &     x101 * y021 *   w1 + x102 * y022 *   w2 + x103 * y023 *   w3
            EFInt(iZeta,iEta,29)=
     &           (y121 *   w1)+       (y122 *   w2)+       (y123 *   w3)
            EFInt(iZeta,iEta,30)=
     &           (y021 * z101)+       (y022 * z102)+       (y023 * z103)
            EFInt(iZeta,iEta,31)=
     &     x201 * y021 *   w1 + x202 * y022 *   w2 + x203 * y023 *   w3
            EFInt(iZeta,iEta,32)=
     &     x101 *(y121 *   w1)+ x102 *(y122 *   w2)+ x103 *(y123 *   w3)
            EFInt(iZeta,iEta,33)=
     &     x101 *(y021 * z101)+ x102 *(y022 * z102)+ x103 *(y023 * z103)
            EFInt(iZeta,iEta,34)=
     &            y221 *   w1 +        y222 *   w2 +        y223 *   w3
            EFInt(iZeta,iEta,35)=
     &            y121 * z101 +        y122 * z102 +        y123 * z103
            EFInt(iZeta,iEta,36)=
     &            y021 * z201 +        y022 * z202 +        y023 * z203
            EFInt(iZeta,iEta,37)=
     &    (x101 * y011)* z011 +(x102 * y012)* z012 +(x103 * y013)* z013
            EFInt(iZeta,iEta,38)=
     &            y111 * z011 +        y112 * z012 +        y113 * z013
            EFInt(iZeta,iEta,39)=
     &            y011 * z111 +        y012 * z112 +        y013 * z113
            EFInt(iZeta,iEta,40)=
     &    (x201 * y011)* z011 +(x202 * y012)* z012 +(x203 * y013)* z013
            EFInt(iZeta,iEta,41)=
     &    (x101 * y111)* z011 +(x102 * y112)* z012 +(x103 * y113)* z013
            EFInt(iZeta,iEta,42)=
     &    (x101 * y011)* z111 +(x102 * y012)* z112 +(x103 * y013)* z113
            EFInt(iZeta,iEta,43)=
     &            y211 * z011 +        y212 * z012 +        y213 * z013
            EFInt(iZeta,iEta,44)=
     &            y111 * z111 +        y112 * z112 +        y113 * z113
            EFInt(iZeta,iEta,45)=
     &            y011 * z211 +        y012 * z212 +        y013 * z213
            EFInt(iZeta,iEta,46)=
     &     x101        * z021 + x102        * z022 + x103        * z023
            EFInt(iZeta,iEta,47)=
     &           (y101 * z021)+       (y102 * z022)+       (y103 * z023)
            EFInt(iZeta,iEta,48)=
     &                   z121 +               z122 +               z123
            EFInt(iZeta,iEta,49)=
     &     x201        * z021 + x202        * z022 + x203        * z023
            EFInt(iZeta,iEta,50)=
     &     x101 *(y101 * z021)+ x102 *(y102 * z022)+ x103 *(y103 * z023)
            EFInt(iZeta,iEta,51)=
     &     x101        * z121 + x102        * z122 + x103        * z123
            EFInt(iZeta,iEta,52)=
     &            y201 * z021 +        y202 * z022 +        y203 * z023
            EFInt(iZeta,iEta,53)=
     &            y101 * z121 +        y102 * z122 +        y103 * z123
            EFInt(iZeta,iEta,54)=
     &                   z221 +               z222 +               z223
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
      PQ2 = PQx**2 + PQy**2 + PQz**2
      Do 14 iEta = 1, nEta
         Do 24 iZeta = 1, nZeta
            ZEInv = One/(Eta(iEta)+Zeta(iZeta)
     >              +(Eta(iEta)*Zeta(iZeta)*ChiI2)*Dble(IsChi))
            rho = Zeta(iZeta)*(Eta(iEta)*ZEInv)
            T = rho * PQ2
            If (T.lt.TMax) Then
               n = iPntr(Int((T+dddx)*xdInv))
               z = T - x0(n)
               w1=(((((CW6(n,1)*z+CW5(n,1))*z+CW4(n,1))*z+CW3(n,1))*z+
     &            CW2(n,1))*z+CW1(n,1))*z+Cw0(n,1)
               w2=(((((CW6(n,2)*z+CW5(n,2))*z+CW4(n,2))*z+CW3(n,2))*z+
     &            CW2(n,2))*z+CW1(n,2))*z+Cw0(n,2)
               w3=(((((CW6(n,3)*z+CW5(n,3))*z+CW4(n,3))*z+CW3(n,3))*z+
     &            CW2(n,3))*z+CW1(n,3))*z+Cw0(n,3)
               r1=(((((CR6(n,1)*z+CR5(n,1))*z+CR4(n,1))*z+CR3(n,1))*z+
     &            CR2(n,1))*z+CR1(n,1))*z+CR0(n,1)
               r2=(((((CR6(n,2)*z+CR5(n,2))*z+CR4(n,2))*z+CR3(n,2))*z+
     &            CR2(n,2))*z+CR1(n,2))*z+CR0(n,2)
               r3=(((((CR6(n,3)*z+CR5(n,3))*z+CR4(n,3))*z+CR3(n,3))*z+
     &            CR2(n,3))*z+CR1(n,3))*z+CR0(n,3)
            Else
               ai = 1.0D0/T
               si = Sqrt(ai)
               w1= HerW(1)*si
               w2= HerW(2)*si
               w3= HerW(3)*si
               r1= HerR2(1)*ai
               r2= HerR2(2)*ai
               r3= HerR2(3)*ai
            End If
            PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv)
            w1 = PreFct*w1
            w2 = PreFct*w2
            w3 = PreFct*w3
            Eu21 =     Eta(iEta)*(r1*ZEInv)
            Eu22 =     Eta(iEta)*(r2*ZEInv)
            Eu23 =     Eta(iEta)*(r3*ZEInv)
            Zu21 =   Zeta(iZeta)*(r1*ZEInv)
            Zu22 =   Zeta(iZeta)*(r2*ZEInv)
            Zu23 =   Zeta(iZeta)*(r3*ZEInv)
            PAQPx1 = - Eu21 * PQx
            PAQPx2 = - Eu22 * PQx
            PAQPx3 = - Eu23 * PQx
            PAQPy1 = - Eu21 * PQy
            PAQPy2 = - Eu22 * PQy
            PAQPy3 = - Eu23 * PQy
            PAQPz1 = - Eu21 * PQz
            PAQPz2 = - Eu22 * PQz
            PAQPz3 = - Eu23 * PQz
            QCPQx1 =  Zu21 * PQx
            QCPQx2 =  Zu22 * PQx
            QCPQx3 =  Zu23 * PQx
            QCPQy1 =  Zu21 * PQy
            QCPQy2 =  Zu22 * PQy
            QCPQy3 =  Zu23 * PQy
            QCPQz1 =  Zu21 * PQz
            QCPQz2 =  Zu22 * PQz
            QCPQz3 =  Zu23 * PQz
            B001 = Half * (r1*ZEInv)
            B002 = Half * (r2*ZEInv)
            B003 = Half * (r3*ZEInv)
            B101 = (Half - Half * Eu21) * ZInv(iZeta)
            B102 = (Half - Half * Eu22) * ZInv(iZeta)
            B103 = (Half - Half * Eu23) * ZInv(iZeta)
            B011 = (Half - Half * Zu21) * EInv(iEta)
            B012 = (Half - Half * Zu22) * EInv(iEta)
            B013 = (Half - Half * Zu23) * EInv(iEta)
            x101= PAQPx1
            x102= PAQPx2
            x103= PAQPx3
            x011= QCPQx1
            x012= QCPQx2
            x013= QCPQx3
            x201= PAQPx1*x101 + B101
            x202= PAQPx2*x102 + B102
            x203= PAQPx3*x103 + B103
            x021= QCPQx1*x011 + B011
            x022= QCPQx2*x012 + B012
            x023= QCPQx3*x013 + B013
            x111= PAQPx1*QCPQx1 + B001
            x112= PAQPx2*QCPQx2 + B002
            x113= PAQPx3*QCPQx3 + B003
            x211= PAQPx1*x111 + B101*x011 + B001*x101
            x212= PAQPx2*x112 + B102*x012 + B002*x102
            x213= PAQPx3*x113 + B103*x013 + B003*x103
            x121= QCPQx1*x111 + B011*x101 + B001*x011
            x122= QCPQx2*x112 + B012*x102 + B002*x012
            x123= QCPQx3*x113 + B013*x103 + B003*x013
            x221= PAQPx1*x121 + B101*x021 + Two*B001*x111
            x222= PAQPx2*x122 + B102*x022 + Two*B002*x112
            x223= PAQPx3*x123 + B103*x023 + Two*B003*x113
            y101= PAQPy1
            y102= PAQPy2
            y103= PAQPy3
            y011= QCPQy1
            y012= QCPQy2
            y013= QCPQy3
            y201= PAQPy1*y101 + B101
            y202= PAQPy2*y102 + B102
            y203= PAQPy3*y103 + B103
            y021= QCPQy1*y011 + B011
            y022= QCPQy2*y012 + B012
            y023= QCPQy3*y013 + B013
            y111= PAQPy1*QCPQy1 + B001
            y112= PAQPy2*QCPQy2 + B002
            y113= PAQPy3*QCPQy3 + B003
            y211= PAQPy1*y111 + B101*y011 + B001*y101
            y212= PAQPy2*y112 + B102*y012 + B002*y102
            y213= PAQPy3*y113 + B103*y013 + B003*y103
            y121= QCPQy1*y111 + B011*y101 + B001*y011
            y122= QCPQy2*y112 + B012*y102 + B002*y012
            y123= QCPQy3*y113 + B013*y103 + B003*y013
            y221= PAQPy1*y121 + B101*y021 + Two*B001*y111
            y222= PAQPy2*y122 + B102*y022 + Two*B002*y112
            y223= PAQPy3*y123 + B103*y023 + Two*B003*y113
            z101= PAQPz1*w1
            z102= PAQPz2*w2
            z103= PAQPz3*w3
            z011= QCPQz1*w1
            z012= QCPQz2*w2
            z013= QCPQz3*w3
            z201= PAQPz1*z101 + B101*w1
            z202= PAQPz2*z102 + B102*w2
            z203= PAQPz3*z103 + B103*w3
            z021= QCPQz1*z011 + B011*w1
            z022= QCPQz2*z012 + B012*w2
            z023= QCPQz3*z013 + B013*w3
            z111= PAQPz1*QCPQz1*w1 + B001*w1
            z112= PAQPz2*QCPQz2*w2 + B002*w2
            z113= PAQPz3*QCPQz3*w3 + B003*w3
            z211= PAQPz1*z111 + B101*z011 + B001*z101
            z212= PAQPz2*z112 + B102*z012 + B002*z102
            z213= PAQPz3*z113 + B103*z013 + B003*z103
            z121= QCPQz1*z111 + B011*z101 + B001*z011
            z122= QCPQz2*z112 + B012*z102 + B002*z012
            z123= QCPQz3*z113 + B013*z103 + B003*z013
            z221= PAQPz1*z121 + B101*z021 + Two*B001*z111
            z222= PAQPz2*z122 + B102*z022 + Two*B002*z112
            z223= PAQPz3*z123 + B103*z023 + Two*B003*z113
            EFInt(iZeta,iEta, 1)=
     &     x221        *   w1 + x222        *   w2 + x223        *   w3
            EFInt(iZeta,iEta, 2)=
     &    (x121 *   w1)* y101 +(x122 *   w2)* y102 +(x123 *   w3)* y103
            EFInt(iZeta,iEta, 3)=
     &     x121        * z101 + x122        * z102 + x123        * z103
            EFInt(iZeta,iEta, 4)=
     &     x021 * y201 *   w1 + x022 * y202 *   w2 + x023 * y203 *   w3
            EFInt(iZeta,iEta, 5)=
     &    (x021 * y101)* z101 +(x022 * y102)* z102 +(x023 * y103)* z103
            EFInt(iZeta,iEta, 6)=
     &     x021        * z201 + x022        * z202 + x023        * z203
            EFInt(iZeta,iEta, 7)=
     &    (x211 *   w1)* y011 +(x212 *   w2)* y012 +(x213 *   w3)* y013
            EFInt(iZeta,iEta, 8)=
     &     x111 * y111 *   w1 + x112 * y112 *   w2 + x113 * y113 *   w3
            EFInt(iZeta,iEta, 9)=
     &     x111 *(y011 * z101)+ x112 *(y012 * z102)+ x113 *(y013 * z103)
            EFInt(iZeta,iEta,10)=
     &     x011 *(y211 *   w1)+ x012 *(y212 *   w2)+ x013 *(y213 *   w3)
            EFInt(iZeta,iEta,11)=
     &    (x011 * y111)* z101 +(x012 * y112)* z102 +(x013 * y113)* z103
            EFInt(iZeta,iEta,12)=
     &    (x011 * y011)* z201 +(x012 * y012)* z202 +(x013 * y013)* z203
            EFInt(iZeta,iEta,13)=
     &     x211        * z011 + x212        * z012 + x213        * z013
            EFInt(iZeta,iEta,14)=
     &     x111 *(y101 * z011)+ x112 *(y102 * z012)+ x113 *(y103 * z013)
            EFInt(iZeta,iEta,15)=
     &     x111        * z111 + x112        * z112 + x113        * z113
            EFInt(iZeta,iEta,16)=
     &    (x011 * y201)* z011 +(x012 * y202)* z012 +(x013 * y203)* z013
            EFInt(iZeta,iEta,17)=
     &    (x011 * y101)* z111 +(x012 * y102)* z112 +(x013 * y103)* z113
            EFInt(iZeta,iEta,18)=
     &     x011        * z211 + x012        * z212 + x013        * z213
            EFInt(iZeta,iEta,19)=
     &     x201 * y021 *   w1 + x202 * y022 *   w2 + x203 * y023 *   w3
            EFInt(iZeta,iEta,20)=
     &     x101 *(y121 *   w1)+ x102 *(y122 *   w2)+ x103 *(y123 *   w3)
            EFInt(iZeta,iEta,21)=
     &     x101 *(y021 * z101)+ x102 *(y022 * z102)+ x103 *(y023 * z103)
            EFInt(iZeta,iEta,22)=
     &            y221 *   w1 +        y222 *   w2 +        y223 *   w3
            EFInt(iZeta,iEta,23)=
     &            y121 * z101 +        y122 * z102 +        y123 * z103
            EFInt(iZeta,iEta,24)=
     &            y021 * z201 +        y022 * z202 +        y023 * z203
            EFInt(iZeta,iEta,25)=
     &    (x201 * y011)* z011 +(x202 * y012)* z012 +(x203 * y013)* z013
            EFInt(iZeta,iEta,26)=
     &    (x101 * y111)* z011 +(x102 * y112)* z012 +(x103 * y113)* z013
            EFInt(iZeta,iEta,27)=
     &    (x101 * y011)* z111 +(x102 * y012)* z112 +(x103 * y013)* z113
            EFInt(iZeta,iEta,28)=
     &            y211 * z011 +        y212 * z012 +        y213 * z013
            EFInt(iZeta,iEta,29)=
     &            y111 * z111 +        y112 * z112 +        y113 * z113
            EFInt(iZeta,iEta,30)=
     &            y011 * z211 +        y012 * z212 +        y013 * z213
            EFInt(iZeta,iEta,31)=
     &     x201        * z021 + x202        * z022 + x203        * z023
            EFInt(iZeta,iEta,32)=
     &     x101 *(y101 * z021)+ x102 *(y102 * z022)+ x103 *(y103 * z023)
            EFInt(iZeta,iEta,33)=
     &     x101        * z121 + x102        * z122 + x103        * z123
            EFInt(iZeta,iEta,34)=
     &            y201 * z021 +        y202 * z022 +        y203 * z023
            EFInt(iZeta,iEta,35)=
     &            y101 * z121 +        y102 * z122 +        y103 * z123
            EFInt(iZeta,iEta,36)=
     &                   z221 +               z222 +               z223
 24      Continue
 14   Continue
*
 99   Continue
*
*     Call qExit('pppp')
      Return
      End
