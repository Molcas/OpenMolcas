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
      Subroutine ppss(EFInt,Zeta,ZInv,nZeta,P,lP,rKappAB,A,B,
     &                       Eta,EInv, nEta,Q,lQ,rKappCD,C,D,
     &                CoorAC,TMax,
     &                iPntr,nPntr,x0,nMax,CW6,CW5,CW4,CW3,CW2,CW1,CW0,
     &                                    CR6,CR5,CR4,CR3,CR2,CR1,CR0,
     &                ddx,HerW,HerR2,IsChi,ChiI2)
************************************************************************
*                                                                      *
* Object: to compute the primitive integrals of type (pp|ss).          *
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
      Real*8 EFInt(nZeta,nEta,9), Zeta(nZeta), Eta(nEta),
     &       CoorAC(3,2), ZInv(nZeta), EInv(nEta),
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
*     Call qEnter('ppss')
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
            Eu21 =   r1*(Eta(iEta)*ZEInv)
            Eu22 =   r2*(Eta(iEta)*ZEInv)
            PAQPx1 = (P(iZeta,1) - CoorAC(1,1)) - Eu21 * PQx
            PAQPx2 = (P(iZeta,1) - CoorAC(1,1)) - Eu22 * PQx
            PAQPy1 = (P(iZeta,2) - CoorAC(2,1)) - Eu21 * PQy
            PAQPy2 = (P(iZeta,2) - CoorAC(2,1)) - Eu22 * PQy
            PAQPz1 = (P(iZeta,3) - CoorAC(3,1)) - Eu21 * PQz
            PAQPz2 = (P(iZeta,3) - CoorAC(3,1)) - Eu22 * PQz
            B101 = (Half - Half * Eu21) * ZInv(iZeta)
            B102 = (Half - Half * Eu22) * ZInv(iZeta)
            x101= PAQPx1
            x102= PAQPx2
            x201= PAQPx1*x101 + B101
            x202= PAQPx2*x102 + B102
            y101= PAQPy1
            y102= PAQPy2
            y201= PAQPy1*y101 + B101
            y202= PAQPy2*y102 + B102
            z101= PAQPz1*w1
            z102= PAQPz2*w2
            z201= PAQPz1*z101 + B101*w1
            z202= PAQPz2*z102 + B102*w2
            PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv)
            EFInt(iZeta,iEta,1) = PreFct * (     x101*w1+     x102*w2)
            EFInt(iZeta,iEta,2) = PreFct * (     y101*w1+     y102*w2)
            EFInt(iZeta,iEta,3) = PreFct * (     z101   +     z102   )
            EFInt(iZeta,iEta,4) = PreFct * (     x201*w1+     x202*w2)
            EFInt(iZeta,iEta,5) = PreFct * (x101*y101*w1+x102*y102*w2)
            EFInt(iZeta,iEta,6) = PreFct * (x101*z101   +x102*z102   )
            EFInt(iZeta,iEta,7) = PreFct * (     y201*w1+     y202*w2)
            EFInt(iZeta,iEta,8) = PreFct * (y101*z101   +y102*z102   )
            EFInt(iZeta,iEta,9) = PreFct * (     z201   +     z202   )
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
            Eu21 =   r1*(Eta(iEta)*ZEInv)
            Eu22 =   r2*(Eta(iEta)*ZEInv)
            B101 = (Half - Half * Eu21) * ZInv(iZeta)
            B102 = (Half - Half * Eu22) * ZInv(iZeta)
            PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv)
            EFInt(iZeta,iEta,1) = (PreFct * (B101*w1+B102*w2))
            EFInt(iZeta,iEta,2) = Zero
            EFInt(iZeta,iEta,3) = Zero
            EFInt(iZeta,iEta,4) = (PreFct * (B101*w1+B102*w2))
            EFInt(iZeta,iEta,5) = Zero
            EFInt(iZeta,iEta,6) = (PreFct * (B101*w1+B102*w2))
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
         PQ2 = (PQx**2 + PQy**2 + PQz**2)
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
            Eu21 =   r1*(Eta(iEta)*ZEInv)
            Eu22 =   r2*(Eta(iEta)*ZEInv)
            PAQPx1 =  - Eu21 * PQx
            PAQPx2 =  - Eu22 * PQx
            PAQPy1 =  - Eu21 * PQy
            PAQPy2 =  - Eu22 * PQy
            PAQPz1 =  - Eu21 * PQz
            PAQPz2 =  - Eu22 * PQz
            B101 = (Half - Half * Eu21) * ZInv(iZeta)
            B102 = (Half - Half * Eu22) * ZInv(iZeta)
            x101= PAQPx1
            x102= PAQPx2
            x201= PAQPx1*x101 + B101
            x202= PAQPx2*x102 + B102
            y101= PAQPy1
            y102= PAQPy2
            y201= PAQPy1*y101 + B101
            y202= PAQPy2*y102 + B102
            z101= PAQPz1*w1
            z102= PAQPz2*w2
            z201= PAQPz1*z101 + B101*w1
            z202= PAQPz2*z102 + B102*w2
            PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv)
            EFInt(iZeta,iEta,1) = PreFct * (     x201*w1+     x202*w2)
            EFInt(iZeta,iEta,2) = PreFct * (x101*y101*w1+x102*y102*w2)
            EFInt(iZeta,iEta,3) = PreFct * (x101*z101   +x102*z102   )
            EFInt(iZeta,iEta,4) = PreFct * (     y201*w1+     y202*w2)
            EFInt(iZeta,iEta,5) = PreFct * (y101*z101   +y102*z102   )
            EFInt(iZeta,iEta,6) = PreFct * (     z201   +     z202   )
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
            Eu21 =   r1*(Eta(iEta)*ZEInv)
            Eu22 =   r2*(Eta(iEta)*ZEInv)
            PAQPx1 = (P(iZeta,1) - CoorAC(1,1)) - Eu21 * PQx
            PAQPx2 = (P(iZeta,1) - CoorAC(1,1)) - Eu22 * PQx
            PAQPy1 = (P(iZeta,2) - CoorAC(2,1)) - Eu21 * PQy
            PAQPy2 = (P(iZeta,2) - CoorAC(2,1)) - Eu22 * PQy
            PAQPz1 = (P(iZeta,3) - CoorAC(3,1)) - Eu21 * PQz
            PAQPz2 = (P(iZeta,3) - CoorAC(3,1)) - Eu22 * PQz
            B101 = (Half - Half * Eu21) * ZInv(iZeta)
            B102 = (Half - Half * Eu22) * ZInv(iZeta)
            x101= PAQPx1
            x102= PAQPx2
            x201= PAQPx1*x101 + B101
            x202= PAQPx2*x102 + B102
            y101= PAQPy1
            y102= PAQPy2
            y201= PAQPy1*y101 + B101
            y202= PAQPy2*y102 + B102
            z101= PAQPz1*w1
            z102= PAQPz2*w2
            z201= PAQPz1*z101 + B101*w1
            z202= PAQPz2*z102 + B102*w2
            PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv)
            EFInt(iZeta,iEta,1) = PreFct * (     x101*w1+     x102*w2)
            EFInt(iZeta,iEta,2) = PreFct * (     y101*w1+     y102*w2)
            EFInt(iZeta,iEta,3) = PreFct * (     z101   +     z102   )
            EFInt(iZeta,iEta,4) = PreFct * (     x201*w1+     x202*w2)
            EFInt(iZeta,iEta,5) = PreFct * (x101*y101*w1+x102*y102*w2)
            EFInt(iZeta,iEta,6) = PreFct * (x101*z101   +x102*z102   )
            EFInt(iZeta,iEta,7) = PreFct * (     y201*w1+     y202*w2)
            EFInt(iZeta,iEta,8) = PreFct * (y101*z101   +y102*z102   )
            EFInt(iZeta,iEta,9) = PreFct * (     z201   +     z202   )
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
            rho = Zeta(iZeta)*(Eta(iEta)*ZEInv)
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
            Eu21 =   r1*(Eta(iEta)*ZEInv)
            Eu22 =   r2*(Eta(iEta)*ZEInv)
            PAQPx1 =  - Eu21 * PQx
            PAQPx2 =  - Eu22 * PQx
            PAQPy1 =  - Eu21 * PQy
            PAQPy2 =  - Eu22 * PQy
            PAQPz1 =  - Eu21 * PQz
            PAQPz2 =  - Eu22 * PQz
            B101 = (Half - Half * Eu21) * ZInv(iZeta)
            B102 = (Half - Half * Eu22) * ZInv(iZeta)
            x101= PAQPx1
            x102= PAQPx2
            x201= PAQPx1*x101 + B101
            x202= PAQPx2*x102 + B102
            y101= PAQPy1
            y102= PAQPy2
            y201= PAQPy1*y101 + B101
            y202= PAQPy2*y102 + B102
            z101= PAQPz1*w1
            z102= PAQPz2*w2
            z201= PAQPz1*z101 + B101*w1
            z202= PAQPz2*z102 + B102*w2
            PreFct = rKappCD(iEta) * rKappAB(iZeta) * Sqrt(ZEInv)
            EFInt(iZeta,iEta,1) = PreFct * (     x201*w1+     x202*w2)
            EFInt(iZeta,iEta,2) = PreFct * (x101*y101*w1+x102*y102*w2)
            EFInt(iZeta,iEta,3) = PreFct * (x101*z101   +x102*z102   )
            EFInt(iZeta,iEta,4) = PreFct * (     y201*w1+     y202*w2)
            EFInt(iZeta,iEta,5) = PreFct * (y101*z101   +y102*z102   )
            EFInt(iZeta,iEta,6) = PreFct * (     z201   +     z202   )
 24      Continue
 14   Continue
*
 99   Continue
*
*     Call qExit('ppss')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(EInv)
      End
