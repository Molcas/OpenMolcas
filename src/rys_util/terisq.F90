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
! Copyright (C) 1990,1992, Roland Lindh                                *
!***********************************************************************

subroutine TERISq(Zeta,Eta,P,Q,rKapab,rKapcd,T,Fact,ZEInv,nT,IsChi,ChiI2,nOrdOp)
!***********************************************************************
!                                                                      *
! Object: to compute entities for the two-electron integrals which are *
!         used in the Rys quadrature to evaluate these integrals.      *
!                                                                      *
!         OBSERVE that the prefactor is only partial!!!                *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             May '92 Modified to fit 2nd order differential scheme for*
!             the gradient estimates.                                  *
!***********************************************************************

use Constants, only: One, Two, Three, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nT, IsChi, nOrdOp
real(kind=wp), intent(in) :: Zeta(nT), Eta(nT), P(nT,3), Q(nT,3), rKapab(nT), rKapcd(nT), ChiI2
real(kind=wp), intent(out) :: T(nT), Fact(nT), ZEInv(nT)
integer(kind=iwp) :: iT
real(kind=wp) :: PQ2, Rho, tmp

#include "macros.fh"
unused_var(Eta)

#ifdef _DEBUGPRINT_
call RecPrt(' Zeta in TERISq',' ',Zeta,nT,1)
call RecPrt(' P in TERISq',' ',P,nT,3)
call RecPrt(' Q in TERISq',' ',Q,nT,3)
call RecPrt(' Kab in TERISq',' ',rKapab,nT,1)
call RecPrt(' Kcd in TERISq',' ',rKapcd,nT,1)
#endif

select case (nOrdOp)

  case (0)

    do iT=1,nT
      tmp = One/(Zeta(iT)+Zeta(iT)+(Zeta(iT)*Zeta(iT)*ChiI2)*real(IsChi,kind=wp))
      ZEInv(iT) = tmp
      Rho = Zeta(iT)*Zeta(iT)*tmp
      PQ2 = (P(iT,1)-Q(iT,1))**2+(P(iT,2)-Q(iT,2))**2+(P(iT,3)-Q(iT,3))**2
      T(iT) = Rho*PQ2
      Fact(iT) = rKapab(iT)*rKapcd(iT)
    end do

  case (1)

    do iT=1,nT
      tmp = One/(Zeta(iT)+Zeta(iT)+(Zeta(iT)*Zeta(iT)*ChiI2)*real(IsChi,kind=wp))
      ZEInv(iT) = tmp
      Rho = Zeta(iT)*Zeta(iT)*tmp
      PQ2 = (P(iT,1)-Q(iT,1))**2+(P(iT,2)-Q(iT,2))**2+(P(iT,3)-Q(iT,3))**2
      T(iT) = Rho*PQ2
      Fact(iT) = rKapab(iT)*rKapcd(iT)*(Two*Rho)
    end do

  case (2)

    do iT=1,nT
      tmp = One/(Zeta(iT)+Zeta(iT)+(Zeta(iT)*Zeta(iT)*ChiI2)*real(IsChi,kind=wp))
      ZEInv(iT) = tmp
      Rho = Zeta(iT)*Zeta(iT)*tmp
      PQ2 = (P(iT,1)-Q(iT,1))**2+(P(iT,2)-Q(iT,2))**2+(P(iT,3)-Q(iT,3))**2
      T(iT) = Rho*PQ2
      Fact(iT) = rKapab(iT)*rKapcd(iT)*(Four*Rho**2/Three)
    end do

end select

#ifdef _DEBUGPRINT_
call RecPrt('Tvalue',' ',T,nT,1)
call RecPrt('Fact  ',' ',Fact,nT,1)
#endif

return

end subroutine TERISq
