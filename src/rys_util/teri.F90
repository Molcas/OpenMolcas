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
! Copyright (C) 1990,2023, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine TERI(Zeta,Eta,P,Q,rKapab,rKapcd,T,Fact,ZEInv,nT,IsChi,ChiI2,nOrdOp)
!***********************************************************************
!                                                                      *
! Object: to entities for the two-electron integrals which are used in *
!         in the Rys quadrature to evaluate these integrals.           *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!***********************************************************************

use Constants, only: One, Two, Three, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nT, IsChi, nOrdOp
real(kind=wp), intent(in) :: Zeta(nT), Eta(nT), P(nT,3), Q(nT,3), rKapab(nT), rKapcd(nT), ChiI2
real(kind=wp), intent(out) :: T(nT), Fact(nT), ZEInv(nT)
integer(kind=iwp) :: iT
real(kind=wp) :: PQ2, Rho, tmp

#ifdef _DEBUGPRINT_
call RecPrt(' Zeta in TERI',' ',Zeta,1,nT)
call RecPrt(' Eta in TERI',' ',Eta,1,nT)
call RecPrt(' P in TERI',' ',P,nT,3)
call RecPrt(' Q in TERI',' ',Q,nT,3)
call RecPrt(' Kab in TERI',' ',rKapab,1,nT)
call RecPrt(' Kcd in TERI',' ',rKapcd,1,nT)
#endif

select case (nOrdOp)

  case (0)

    do iT=1,nT
      tmp = One/(Zeta(iT)+Eta(iT)+(Eta(iT)*Zeta(iT)*ChiI2)*real(IsChi,kind=wp))
      ZEInv(iT) = tmp
      Rho = Zeta(iT)*Eta(iT)*tmp
      PQ2 = (P(iT,1)-Q(iT,1))**2+(P(iT,2)-Q(iT,2))**2+(P(iT,3)-Q(iT,3))**2
      T(iT) = Rho*PQ2
      Fact(iT) = rKapab(iT)*rKapcd(iT)*sqrt(tmp)
    end do

  case (1)

    do iT=1,nT
      tmp = One/(Zeta(iT)+Eta(iT)+(Eta(iT)*Zeta(iT)*ChiI2)*real(IsChi,kind=wp))
      ZEInv(iT) = tmp
      Rho = Zeta(iT)*Eta(iT)*tmp
      PQ2 = (P(iT,1)-Q(iT,1))**2+(P(iT,2)-Q(iT,2))**2+(P(iT,3)-Q(iT,3))**2
      T(iT) = Rho*PQ2
      Fact(iT) = rKapab(iT)*rKapcd(iT)*sqrt(tmp)*Two*Rho
    end do

  case (2)

    do iT=1,nT
      tmp = One/(Zeta(iT)+Eta(iT)+(Eta(iT)*Zeta(iT)*ChiI2)*real(IsChi,kind=wp))
      ZEInv(iT) = tmp
      Rho = Zeta(iT)*Eta(iT)*tmp
      PQ2 = (P(iT,1)-Q(iT,1))**2+(P(iT,2)-Q(iT,2))**2+(P(iT,3)-Q(iT,3))**2
      T(iT) = Rho*PQ2
      Fact(iT) = rKapab(iT)*rKapcd(iT)*sqrt(tmp)*(Four*Rho**2/Three)
    end do

end select
#ifdef _DEBUGPRINT_
call RecPrt('Tvalue',' ',T,1,nT)
call RecPrt('Fact  ',' ',Fact,1,nT)
#endif

return

end subroutine TERI
