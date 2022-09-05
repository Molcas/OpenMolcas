!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994, Roland Lindh                                     *
!***********************************************************************

subroutine ppss(EFInt,Zeta,ZInv,nZeta,P,lP,rKappAB,A,B,Eta,EInv,nEta,Q,lQ,rKappCD,C,D,CoorAC,TMax,iPntr,nPntr,x0,nMax,CW6,CW5,CW4, &
                CW3,CW2,CW1,CW0,CR6,CR5,CR4,CR3,CR2,CR1,CR0,ddx,HerW,HerR2,IsChi,ChiI2)
!***********************************************************************
!                                                                      *
! Object: to compute the primitive integrals of type (pp|ss).          *
!                                                                      *
!  Author:    Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN. 1994                                    *
!***********************************************************************

use Constants, only: Zero, One, Ten, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, lP, nEta, lQ, nPntr, iPntr(nPntr), nMax, IsChi
real(kind=wp), intent(out) :: EFInt(nZeta,nEta,9)
real(kind=wp), intent(in) :: Zeta(nZeta), ZInv(nZeta), P(lP,3), rKappAB(nZeta), A(3), B(3), Eta(nEta), EInv(nEta), Q(lQ,3), &
                             rKappCD(nEta), C(3), D(3), CoorAC(3,2), TMax, x0(nMax), CW6(nMax,2), CW5(nMax,2), CW4(nMax,2), &
                             CW3(nMax,2), CW2(nMax,2), CW1(nMax,2), CW0(nMax,2), CR6(nMax,2), CR5(nMax,2), CR4(nMax,2), &
                             CR3(nMax,2), CR2(nMax,2), CR1(nMax,2), CR0(nMax,2), ddx, HerW(2), HerR2(2), ChiI2
integer(kind=iwp) :: iEta, iZeta, n
real(kind=wp) :: ai, B101, B102, dddx, Eu21, Eu22, PAQPx1, PAQPx2, PAQPy1, PAQPy2, PAQPz1, PAQPz2, PQ2, PQx, PQy, PQz, PreFct, r1, &
                 r2, rho, si, T, w1, w2, x101, x102, x201, x202, xdInv, y101, y102, y201, y202, z, z101, z102, z201, z202, ZEInv
logical(kind=iwp) :: ABeqCD, EQ

#include "macros.fh"
unused_var(EInv)

xdInv = One/ddx
dddx = ddx/Ten+ddx

ABeqCD = EQ(A,B) .and. EQ(A,C) .and. EQ(A,D)

if (ABeqCD) then

  ! AAAA case

  z = -x0(1)
  w1 = (((((CW6(1,1)*z+CW5(1,1))*z+CW4(1,1))*z+CW3(1,1))*z+CW2(1,1))*z+CW1(1,1))*z+Cw0(1,1)
  w2 = (((((CW6(1,2)*z+CW5(1,2))*z+CW4(1,2))*z+CW3(1,2))*z+CW2(1,2))*z+CW1(1,2))*z+Cw0(1,2)
  r1 = (((((CR6(1,1)*z+CR5(1,1))*z+CR4(1,1))*z+CR3(1,1))*z+CR2(1,1))*z+CR1(1,1))*z+CR0(1,1)
  r2 = (((((CR6(1,2)*z+CR5(1,2))*z+CR4(1,2))*z+CR3(1,2))*z+CR2(1,2))*z+CR1(1,2))*z+CR0(1,2)
  do iEta=1,nEta
    do iZeta=1,nZeta
      ZEInv = One/(Eta(iEta)+Zeta(iZeta)+(Eta(iEta)*Zeta(iZeta)*ChiI2)*real(IsChi,kind=wp))
      Eu21 = r1*(Eta(iEta)*ZEInv)
      Eu22 = r2*(Eta(iEta)*ZEInv)
      B101 = (Half-Half*Eu21)*ZInv(iZeta)
      B102 = (Half-Half*Eu22)*ZInv(iZeta)
      PreFct = rKappCD(iEta)*rKappAB(iZeta)*sqrt(ZEInv)
      EFInt(iZeta,iEta,1) = (PreFct*(B101*w1+B102*w2))
      EFInt(iZeta,iEta,2) = Zero
      EFInt(iZeta,iEta,3) = Zero
      EFInt(iZeta,iEta,4) = (PreFct*(B101*w1+B102*w2))
      EFInt(iZeta,iEta,5) = Zero
      EFInt(iZeta,iEta,6) = (PreFct*(B101*w1+B102*w2))
    end do
  end do

else if (EQ(A,B) .and. (.not. EQ(C,D))) then

  ! AACD case

  do iEta=1,nEta
    do iZeta=1,nZeta
      PQx = CoorAC(1,1)-Q(iEta,1)
      PQy = CoorAC(2,1)-Q(iEta,2)
      PQz = CoorAC(3,1)-Q(iEta,3)
      PQ2 = (PQx**2+PQy**2+PQz**2)
      ZEInv = One/(Eta(iEta)+Zeta(iZeta)+(Eta(iEta)*Zeta(iZeta)*ChiI2)*real(IsChi,kind=wp))
      rho = Zeta(iZeta)*(Eta(iEta)*ZEInv)
      T = rho*PQ2
      if (T < TMax) then
        n = iPntr(int((T+dddx)*xdInv))
        z = T-x0(n)
        w1 = (((((CW6(n,1)*z+CW5(n,1))*z+CW4(n,1))*z+CW3(n,1))*z+CW2(n,1))*z+CW1(n,1))*z+Cw0(n,1)
        w2 = (((((CW6(n,2)*z+CW5(n,2))*z+CW4(n,2))*z+CW3(n,2))*z+CW2(n,2))*z+CW1(n,2))*z+Cw0(n,2)
        r1 = (((((CR6(n,1)*z+CR5(n,1))*z+CR4(n,1))*z+CR3(n,1))*z+CR2(n,1))*z+CR1(n,1))*z+CR0(n,1)
        r2 = (((((CR6(n,2)*z+CR5(n,2))*z+CR4(n,2))*z+CR3(n,2))*z+CR2(n,2))*z+CR1(n,2))*z+CR0(n,2)
      else
        ai = One/T
        si = sqrt(ai)
        w1 = HerW(1)*si
        w2 = HerW(2)*si
        r1 = HerR2(1)*ai
        r2 = HerR2(2)*ai
      end if
      Eu21 = r1*(Eta(iEta)*ZEInv)
      Eu22 = r2*(Eta(iEta)*ZEInv)
      PAQPx1 = -Eu21*PQx
      PAQPx2 = -Eu22*PQx
      PAQPy1 = -Eu21*PQy
      PAQPy2 = -Eu22*PQy
      PAQPz1 = -Eu21*PQz
      PAQPz2 = -Eu22*PQz
      B101 = (Half-Half*Eu21)*ZInv(iZeta)
      B102 = (Half-Half*Eu22)*ZInv(iZeta)
      x101 = PAQPx1
      x102 = PAQPx2
      x201 = PAQPx1*x101+B101
      x202 = PAQPx2*x102+B102
      y101 = PAQPy1
      y102 = PAQPy2
      y201 = PAQPy1*y101+B101
      y202 = PAQPy2*y102+B102
      z101 = PAQPz1*w1
      z102 = PAQPz2*w2
      z201 = PAQPz1*z101+B101*w1
      z202 = PAQPz2*z102+B102*w2
      PreFct = rKappCD(iEta)*rKappAB(iZeta)*sqrt(ZEInv)
      EFInt(iZeta,iEta,1) = PreFct*(x201*w1+x202*w2)
      EFInt(iZeta,iEta,2) = PreFct*(x101*y101*w1+x102*y102*w2)
      EFInt(iZeta,iEta,3) = PreFct*(x101*z101+x102*z102)
      EFInt(iZeta,iEta,4) = PreFct*(y201*w1+y202*w2)
      EFInt(iZeta,iEta,5) = PreFct*(y101*z101+y102*z102)
      EFInt(iZeta,iEta,6) = PreFct*(z201+z202)
    end do
  end do

else if ((.not. EQ(A,B)) .and. EQ(C,D)) then

  ! ABCC case

  do iEta=1,nEta
    do iZeta=1,nZeta
      ZEInv = One/(Eta(iEta)+Zeta(iZeta)+(Eta(iEta)*Zeta(iZeta)*ChiI2)*real(IsChi,kind=wp))
      rho = Zeta(iZeta)*(Eta(iEta)*ZEInv)
      PQx = P(iZeta,1)-CoorAC(1,2)
      PQy = P(iZeta,2)-CoorAC(2,2)
      PQz = P(iZeta,3)-CoorAC(3,2)
      T = rho*(PQx**2+PQy**2+PQz**2)
      if (T < TMax) then
        n = iPntr(int((T+dddx)*xdInv))
        z = T-x0(n)
        w1 = (((((CW6(n,1)*z+CW5(n,1))*z+CW4(n,1))*z+CW3(n,1))*z+CW2(n,1))*z+CW1(n,1))*z+Cw0(n,1)
        w2 = (((((CW6(n,2)*z+CW5(n,2))*z+CW4(n,2))*z+CW3(n,2))*z+CW2(n,2))*z+CW1(n,2))*z+Cw0(n,2)
        r1 = (((((CR6(n,1)*z+CR5(n,1))*z+CR4(n,1))*z+CR3(n,1))*z+CR2(n,1))*z+CR1(n,1))*z+CR0(n,1)
        r2 = (((((CR6(n,2)*z+CR5(n,2))*z+CR4(n,2))*z+CR3(n,2))*z+CR2(n,2))*z+CR1(n,2))*z+CR0(n,2)
      else
        ai = One/T
        si = sqrt(ai)
        w1 = HerW(1)*si
        w2 = HerW(2)*si
        r1 = HerR2(1)*ai
        r2 = HerR2(2)*ai
      end if
      Eu21 = r1*(Eta(iEta)*ZEInv)
      Eu22 = r2*(Eta(iEta)*ZEInv)
      PAQPx1 = (P(iZeta,1)-CoorAC(1,1))-Eu21*PQx
      PAQPx2 = (P(iZeta,1)-CoorAC(1,1))-Eu22*PQx
      PAQPy1 = (P(iZeta,2)-CoorAC(2,1))-Eu21*PQy
      PAQPy2 = (P(iZeta,2)-CoorAC(2,1))-Eu22*PQy
      PAQPz1 = (P(iZeta,3)-CoorAC(3,1))-Eu21*PQz
      PAQPz2 = (P(iZeta,3)-CoorAC(3,1))-Eu22*PQz
      B101 = (Half-Half*Eu21)*ZInv(iZeta)
      B102 = (Half-Half*Eu22)*ZInv(iZeta)
      x101 = PAQPx1
      x102 = PAQPx2
      x201 = PAQPx1*x101+B101
      x202 = PAQPx2*x102+B102
      y101 = PAQPy1
      y102 = PAQPy2
      y201 = PAQPy1*y101+B101
      y202 = PAQPy2*y102+B102
      z101 = PAQPz1*w1
      z102 = PAQPz2*w2
      z201 = PAQPz1*z101+B101*w1
      z202 = PAQPz2*z102+B102*w2
      PreFct = rKappCD(iEta)*rKappAB(iZeta)*sqrt(ZEInv)
      EFInt(iZeta,iEta,1) = PreFct*(x101*w1+x102*w2)
      EFInt(iZeta,iEta,2) = PreFct*(y101*w1+y102*w2)
      EFInt(iZeta,iEta,3) = PreFct*(z101+z102)
      EFInt(iZeta,iEta,4) = PreFct*(x201*w1+x202*w2)
      EFInt(iZeta,iEta,5) = PreFct*(x101*y101*w1+x102*y102*w2)
      EFInt(iZeta,iEta,6) = PreFct*(x101*z101+x102*z102)
      EFInt(iZeta,iEta,7) = PreFct*(y201*w1+y202*w2)
      EFInt(iZeta,iEta,8) = PreFct*(y101*z101+y102*z102)
      EFInt(iZeta,iEta,9) = PreFct*(z201+z202)
    end do
  end do

else if (EQ(A,B) .and. EQ(C,D)) then

  ! AACC case

  PQx = CoorAC(1,1)-CoorAC(1,2)
  PQy = CoorAC(2,1)-CoorAC(2,2)
  PQz = CoorAC(3,1)-CoorAC(3,2)
  PQ2 = (PQx**2+PQy**2+PQz**2)
  do iEta=1,nEta
    do iZeta=1,nZeta
      ZEInv = One/(Eta(iEta)+Zeta(iZeta)+(Eta(iEta)*Zeta(iZeta)*ChiI2)*real(IsChi,kind=wp))
      rho = Zeta(iZeta)*(Eta(iEta)*ZEInv)
      T = rho*PQ2
      if (T < TMax) then
        n = iPntr(int((T+dddx)*xdInv))
        z = T-x0(n)
        w1 = (((((CW6(n,1)*z+CW5(n,1))*z+CW4(n,1))*z+CW3(n,1))*z+CW2(n,1))*z+CW1(n,1))*z+Cw0(n,1)
        w2 = (((((CW6(n,2)*z+CW5(n,2))*z+CW4(n,2))*z+CW3(n,2))*z+CW2(n,2))*z+CW1(n,2))*z+Cw0(n,2)
        r1 = (((((CR6(n,1)*z+CR5(n,1))*z+CR4(n,1))*z+CR3(n,1))*z+CR2(n,1))*z+CR1(n,1))*z+CR0(n,1)
        r2 = (((((CR6(n,2)*z+CR5(n,2))*z+CR4(n,2))*z+CR3(n,2))*z+CR2(n,2))*z+CR1(n,2))*z+CR0(n,2)
      else
        ai = One/T
        si = sqrt(ai)
        w1 = HerW(1)*si
        w2 = HerW(2)*si
        r1 = HerR2(1)*ai
        r2 = HerR2(2)*ai
      end if
      Eu21 = r1*(Eta(iEta)*ZEInv)
      Eu22 = r2*(Eta(iEta)*ZEInv)
      PAQPx1 = -Eu21*PQx
      PAQPx2 = -Eu22*PQx
      PAQPy1 = -Eu21*PQy
      PAQPy2 = -Eu22*PQy
      PAQPz1 = -Eu21*PQz
      PAQPz2 = -Eu22*PQz
      B101 = (Half-Half*Eu21)*ZInv(iZeta)
      B102 = (Half-Half*Eu22)*ZInv(iZeta)
      x101 = PAQPx1
      x102 = PAQPx2
      x201 = PAQPx1*x101+B101
      x202 = PAQPx2*x102+B102
      y101 = PAQPy1
      y102 = PAQPy2
      y201 = PAQPy1*y101+B101
      y202 = PAQPy2*y102+B102
      z101 = PAQPz1*w1
      z102 = PAQPz2*w2
      z201 = PAQPz1*z101+B101*w1
      z202 = PAQPz2*z102+B102*w2
      PreFct = rKappCD(iEta)*rKappAB(iZeta)*sqrt(ZEInv)
      EFInt(iZeta,iEta,1) = PreFct*(x201*w1+x202*w2)
      EFInt(iZeta,iEta,2) = PreFct*(x101*y101*w1+x102*y102*w2)
      EFInt(iZeta,iEta,3) = PreFct*(x101*z101+x102*z102)
      EFInt(iZeta,iEta,4) = PreFct*(y201*w1+y202*w2)
      EFInt(iZeta,iEta,5) = PreFct*(y101*z101+y102*z102)
      EFInt(iZeta,iEta,6) = PreFct*(z201+z202)
    end do
  end do

else

  ! ABCD case

  do iEta=1,nEta
    do iZeta=1,nZeta
      ZEInv = One/(Eta(iEta)+Zeta(iZeta)+(Eta(iEta)*Zeta(iZeta)*ChiI2)*real(IsChi,kind=wp))
      rho = Zeta(iZeta)*(Eta(iEta)*ZEInv)
      PQx = P(iZeta,1)-Q(iEta,1)
      PQy = P(iZeta,2)-Q(iEta,2)
      PQz = P(iZeta,3)-Q(iEta,3)
      T = rho*(PQx**2+PQy**2+PQz**2)
      if (T < TMax) then
        n = iPntr(int((T+dddx)*xdInv))
        z = T-x0(n)
        w1 = (((((CW6(n,1)*z+CW5(n,1))*z+CW4(n,1))*z+CW3(n,1))*z+CW2(n,1))*z+CW1(n,1))*z+Cw0(n,1)
        w2 = (((((CW6(n,2)*z+CW5(n,2))*z+CW4(n,2))*z+CW3(n,2))*z+CW2(n,2))*z+CW1(n,2))*z+Cw0(n,2)
        r1 = (((((CR6(n,1)*z+CR5(n,1))*z+CR4(n,1))*z+CR3(n,1))*z+CR2(n,1))*z+CR1(n,1))*z+CR0(n,1)
        r2 = (((((CR6(n,2)*z+CR5(n,2))*z+CR4(n,2))*z+CR3(n,2))*z+CR2(n,2))*z+CR1(n,2))*z+CR0(n,2)
      else
        ai = One/T
        si = sqrt(ai)
        w1 = HerW(1)*si
        w2 = HerW(2)*si
        r1 = HerR2(1)*ai
        r2 = HerR2(2)*ai
      end if
      Eu21 = r1*(Eta(iEta)*ZEInv)
      Eu22 = r2*(Eta(iEta)*ZEInv)
      PAQPx1 = (P(iZeta,1)-CoorAC(1,1))-Eu21*PQx
      PAQPx2 = (P(iZeta,1)-CoorAC(1,1))-Eu22*PQx
      PAQPy1 = (P(iZeta,2)-CoorAC(2,1))-Eu21*PQy
      PAQPy2 = (P(iZeta,2)-CoorAC(2,1))-Eu22*PQy
      PAQPz1 = (P(iZeta,3)-CoorAC(3,1))-Eu21*PQz
      PAQPz2 = (P(iZeta,3)-CoorAC(3,1))-Eu22*PQz
      B101 = (Half-Half*Eu21)*ZInv(iZeta)
      B102 = (Half-Half*Eu22)*ZInv(iZeta)
      x101 = PAQPx1
      x102 = PAQPx2
      x201 = PAQPx1*x101+B101
      x202 = PAQPx2*x102+B102
      y101 = PAQPy1
      y102 = PAQPy2
      y201 = PAQPy1*y101+B101
      y202 = PAQPy2*y102+B102
      z101 = PAQPz1*w1
      z102 = PAQPz2*w2
      z201 = PAQPz1*z101+B101*w1
      z202 = PAQPz2*z102+B102*w2
      PreFct = rKappCD(iEta)*rKappAB(iZeta)*sqrt(ZEInv)
      EFInt(iZeta,iEta,1) = PreFct*(x101*w1+x102*w2)
      EFInt(iZeta,iEta,2) = PreFct*(y101*w1+y102*w2)
      EFInt(iZeta,iEta,3) = PreFct*(z101+z102)
      EFInt(iZeta,iEta,4) = PreFct*(x201*w1+x202*w2)
      EFInt(iZeta,iEta,5) = PreFct*(x101*y101*w1+x102*y102*w2)
      EFInt(iZeta,iEta,6) = PreFct*(x101*z101+x102*z102)
      EFInt(iZeta,iEta,7) = PreFct*(y201*w1+y202*w2)
      EFInt(iZeta,iEta,8) = PreFct*(y101*z101+y102*z102)
      EFInt(iZeta,iEta,9) = PreFct*(z201+z202)
    end do
  end do

end if

return

end subroutine ppss
