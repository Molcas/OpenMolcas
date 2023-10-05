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
! Copyright (C) 1990,1991, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine TERIS(Zeta,Eta,P,Q,rKapab,rKapcd,T,Fact,ZEInv,nT,IsChi,ChiI2,nOrdOp)
!***********************************************************************
!                                                                      *
! Object: compute the arguments for the reduced list of integrals which*
!         are used in prescreening.                                    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             June '91, modified for k2 loop.                          *
!***********************************************************************

use Constants, only: Zero, One, Two, Three, Four, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nT, IsChi, nOrdOp
real(kind=wp), intent(in) :: Zeta(nT), Eta(nT), P(nT,3), Q(nT,3), rKapab(nT), rKapcd(nT), ChiI2
real(kind=wp), intent(out) :: T(nT), Fact(nT), ZEInv(nT)
integer(kind=iwp) :: iT
real(kind=wp) :: Rho, tmp

#include "macros.fh"
unused_var(Eta)
unused_var(P)
unused_var(Q)
unused_var(rKapcd)

#ifdef _DEBUGPRINT_
call RecPrt(' Zeta in TERIS',' ',Zeta,nT,1)
call RecPrt(' Kab in TERIS',' ',rKapab,nT,1)
#endif

T(:) = Zero

select case (nOrdOp)

  case (0)

    do iT=1,nT
      tmp = One/(Zeta(iT)+Zeta(iT)+(Zeta(iT)*Zeta(iT)*ChiI2)*real(IsChi,kind=wp))
      ZEInv(iT) = tmp
      Fact(iT) = rKapab(iT)**2*sqrt(tmp)
    end do

  case (1)

    do iT=1,nT
      tmp = One/(Zeta(iT)+Zeta(iT)+(Zeta(iT)*Zeta(iT)*ChiI2)*real(IsChi,kind=wp))
      ZEInv(iT) = tmp
      Rho = Zeta(iT)*Half
      Fact(iT) = rKapab(iT)**2*sqrt(tmp)*Two*Rho
    end do

  case (2)

    do iT=1,nT
      tmp = One/(Zeta(iT)+Zeta(iT)+(Zeta(iT)*Zeta(iT)*ChiI2)*real(IsChi,kind=wp))
      ZEInv(iT) = tmp
      Rho = Zeta(iT)*Half
      Fact(iT) = rKapab(iT)**2*sqrt(tmp)*(Four*Rho**2/Three)
    end do

end select
#ifdef _DEBUGPRINT_
call RecPrt('In TERIS: Tvalue',' ',T,nT,1)
call RecPrt('In TERIS: Fact  ',' ',Fact,nT,1)
#endif

return

end subroutine TERIS
