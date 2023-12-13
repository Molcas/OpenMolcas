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

subroutine ssss(EFInt,Zeta,nZeta,P,lP,rKappAB,A,B,Eta,nEta,Q,lQ,rKappCD,C,D,TMax,iPntr,nPntr,x0,nMax,W6,W5,W4,W3,W2,W1,W0,ddx, &
                HerW,IsChi,ChiI2)
!***********************************************************************
!                                                                      *
! Object: to compute the primitive integrals of type (ss|ss).          *
!                                                                      *
!  Author:    Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN. 1994                                    *
!***********************************************************************

use Constants, only: One, Ten
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, lP, nEta, lQ, nPntr, iPntr(nPntr), nMax, IsChi
real(kind=wp), intent(out) :: EFInt(nZeta,nEta)
real(kind=wp), intent(in) :: Zeta(nZeta), P(lP,3), rKappAB(nZeta), A(3), B(3), Eta(nEta), Q(lQ,3), rKappCD(nEta), C(3), D(3), &
                             TMax, x0(nMax), W6(nMax), W5(nMax), W4(nMax), W3(nMax), W2(nMax), W1(nMax), W0(nMax), ddx, HerW, ChiI2
integer(kind=iwp) :: iEta, iZeta, n
real(kind=wp) :: dddx, PQ2, rho, T, w, xdInv, z, ZE, ZEInv
logical(kind=iwp) :: ABeqCD
logical(kind=iwp), external :: EQ

xdInv = One/ddx
dddx = ddx/Ten+ddx

ABeqCD = EQ(A,B) .and. EQ(A,C) .and. EQ(A,D)

if (ABeqCD) then

  z = -x0(1)
  w = (((((W6(1)*z+W5(1))*z+W4(1))*z+W3(1))*z+W2(1))*z+W1(1))*z+w0(1)
  do iEta=1,nEta
    do iZeta=1,nZeta
      ZEInv = One/(Eta(iEta)+Zeta(iZeta)+(Eta(iEta)*Zeta(iZeta)*ChiI2)*real(IsChi,kind=wp))
      EFInt(iZeta,iEta) = rKappCD(iEta)*rKappAB(iZeta)*sqrt(ZEInv)*w
    end do
  end do

else

  do iEta=1,nEta
    do iZeta=1,nZeta
      ZEInv = One/(Eta(iEta)+Zeta(iZeta)+(Eta(iEta)*Zeta(iZeta)*ChiI2)*real(IsChi,kind=wp))
      ZE = Zeta(iZeta)*Eta(iEta)
      rho = ZE*ZEInv
      PQ2 = (P(iZeta,1)-Q(iEta,1))**2+(P(iZeta,2)-Q(iEta,2))**2+(P(iZeta,3)-Q(iEta,3))**2
      T = rho*PQ2
      if (T < TMax) then
        n = iPntr(int((T+dddx)*xdInv))
        z = T-x0(n)
        w = (((((W6(n)*z+W5(n))*z+W4(n))*z+W3(n))*z+W2(n))*z+W1(n))*z+w0(n)
        EFInt(iZeta,iEta) = rKappCD(iEta)*rKappAB(iZeta)*sqrt(ZEInv)*w
      else
        EFInt(iZeta,iEta) = rKappCD(iEta)*rKappAB(iZeta)*HerW*sqrt(One/(ZE*PQ2))
      end if
    end do
  end do

end if

return

end subroutine ssss
