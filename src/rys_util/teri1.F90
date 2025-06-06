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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine TERI1( &
#                define _CALLING_
#                include "tval1_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: to entities for the two-electron integrals which are used in *
!         in the Rys quadrature to evaluate these integrals.           *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!***********************************************************************

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
#include "tval1_interface.fh"
integer(kind=iwp) :: iT
real(kind=wp) :: PQ2, Rho, tmp

#ifdef _DEBUGPRINT_
call RecPrt(' Zeta in TERI1',' ',Zeta,nT,1)
call RecPrt(' Eta in TERI1',' ',Eta,nT,1)
call RecPrt(' P in TERI1',' ',P,nT,3)
call RecPrt(' Q in TERI1',' ',Q,nT,3)
#endif

do iT=1,nT
  tmp = One/(Zeta(iT)+Eta(iT)+(Eta(iT)*Zeta(iT)*ChiI2)*real(IsChi,kind=wp))
  ZEInv(iT) = tmp
  Rho = Zeta(iT)*Eta(iT)*tmp
  PQ2 = (P(iT,1)-Q(iT,1))**2+(P(iT,2)-Q(iT,2))**2+(P(iT,3)-Q(iT,3))**2
  T(iT) = Rho*PQ2
end do

#ifdef _DEBUGPRINT_
call RecPrt('Tvalue',' ',T,nT,1)
#endif

return

end subroutine TERI1
