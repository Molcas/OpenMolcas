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

subroutine psss(EFInt,Zeta,nZeta,P,lP,rKappAB,A,B,Eta,nEta,Q,lQ,rKappCD,C,D,CoorAC,TMax,iPntr,nPntr,x0,nMax,W6,W5,W4,W3,W2,W1,W0, &
                R6,R5,R4,R3,R2,R1,R0,ddx,HerW,HerR2,IsChi,ChiI2)
!***********************************************************************
!                                                                      *
! Object: to compute the primitive integrals of type (ps|ss).          *
!                                                                      *
!  Author:    Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN. 1994                                    *
!***********************************************************************

use Constants, only: Zero, One, Ten
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, lP, nEta, lQ, nPntr, iPntr(nPntr), nMax, IsChi
real(kind=wp), intent(out) :: EFInt(nZeta,nEta,3)
real(kind=wp), intent(in) :: Zeta(nZeta), P(lP,3), rKappAB(nZeta), A(3), B(3), Eta(nEta), Q(lQ,3), rKappCD(nEta), C(3), D(3), &
                             CoorAC(3,2), TMax, x0(nMax), W6(nMax), W5(nMax), W4(nMax), W3(nMax), W2(nMax), W1(nMax), W0(nMax), &
                             R6(nMax), R5(nMax), R4(nMax), R3(nMax), R2(nMax), R1(nMax), R0(nMax), ddx, HerW, HerR2, ChiI2
integer(kind=iwp) :: iEta, iZeta, n
real(kind=wp) :: dddx, Eu2, PAQPx, PAQPy, PAQPz, PQ2, PQx, PQy, PQz, PreFct, r, rho, t, w, xdInv, z, ZE, ZEInv
logical(kind=iwp) :: ABeqCD, EQ

xdInv = One/ddx
dddx = ddx/Ten+ddx

ABeqCD = EQ(A,B) .and. EQ(A,C) .and. EQ(A,D)

if (ABeqCD) then

  ! CCCC case

  EFInt(:,:,:) = Zero

else if (EQ(A,B)) then

  ! AACD case

  do iEta=1,nEta
    do iZeta=1,nZeta
      PQx = (Q(iEta,1)-CoorAC(1,1))
      PQy = (Q(iEta,2)-CoorAC(2,1))
      PQz = (Q(iEta,3)-CoorAC(3,1))
      PQ2 = PQx**2+PQy**2+PQz**2
      ZEInv = One/(Eta(iEta)+Zeta(iZeta)+(Eta(iEta)*Zeta(iZeta)*ChiI2)*real(IsChi,kind=wp))
      ZE = Zeta(iZeta)*Eta(iEta)
      rho = ZE*ZEInv
      T = rho*PQ2
      if (T < TMax) then
        n = iPntr(int((T+dddx)*xdInv))
        z = T-x0(n)
        w = (((((W6(n)*z+W5(n))*z+W4(n))*z+W3(n))*z+W2(n))*z+W1(n))*z+w0(n)
        r = (((((R6(n)*z+R5(n))*z+R4(n))*z+R3(n))*z+R2(n))*z+R1(n))*z+R0(n)
        Eu2 = r*(Eta(iEta)*ZEInv)
        PreFct = rKappCD(iEta)*rKappAB(iZeta)*sqrt(ZEInv)*w
      else
        Eu2 = HerR2/(Zeta(iZeta)*PQ2)
        PreFct = rKappCD(iEta)*rKappAB(iZeta)*HerW/sqrt(ZE*PQ2)
      end if
      EFInt(iZeta,iEta,1) = PreFct*Eu2*PQx
      EFInt(iZeta,iEta,2) = PreFct*Eu2*PQy
      EFInt(iZeta,iEta,3) = PreFct*Eu2*PQz
    end do
  end do

else

  ! ABCD case

  do iEta=1,nEta
    do iZeta=1,nZeta
      ZEInv = One/(Eta(iEta)+Zeta(iZeta)+(Eta(iEta)*Zeta(iZeta)*ChiI2)*real(IsChi,kind=wp))
      ZE = Zeta(iZeta)*Eta(iEta)
      rho = ZE*ZEInv
      PQx = P(iZeta,1)-Q(iEta,1)
      PQy = P(iZeta,2)-Q(iEta,2)
      PQz = P(iZeta,3)-Q(iEta,3)
      PQ2 = PQx**2+PQy**2+PQz**2
      T = rho*PQ2
      if (T < TMax) then
        n = iPntr(int((T+dddx)*xdInv))
        z = T-x0(n)
        w = (((((W6(n)*z+W5(n))*z+W4(n))*z+W3(n))*z+W2(n))*z+W1(n))*z+w0(n)
        r = (((((R6(n)*z+R5(n))*z+R4(n))*z+R3(n))*z+R2(n))*z+R1(n))*z+R0(n)
        Eu2 = r*(Eta(iEta)*ZEInv)
        PreFct = rKappCD(iEta)*rKappAB(iZeta)*sqrt(ZEInv)*w
      else
        Eu2 = HerR2/(Zeta(iZeta)*PQ2)
        PreFct = rKappCD(iEta)*rKappAB(iZeta)*HerW/sqrt(ZE*PQ2)
      end if
      PAQPx = P(iZeta,1)-CoorAC(1,1)-Eu2*PQx
      PAQPy = P(iZeta,2)-CoorAC(2,1)-Eu2*PQy
      PAQPz = P(iZeta,3)-CoorAC(3,1)-Eu2*PQz
      EFInt(iZeta,iEta,1) = PreFct*PAQPx
      EFInt(iZeta,iEta,2) = PreFct*PAQPy
      EFInt(iZeta,iEta,3) = PreFct*PAQPz
    end do
  end do

end if

return

end subroutine psss
