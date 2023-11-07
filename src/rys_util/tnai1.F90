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
!               1990, IBM                                              *
!***********************************************************************

subroutine TNAI1(Zeta,Eta,P,Q,nT,T,ZEInv,IsChi,ChiI2)
!***********************************************************************
!                                                                      *
! Object: to compute entities for the nuclear attraction integrals     *
!         which are used in the Rys quadrature to evaluate these       *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             March '92 modified to gradient calculation.              *
!***********************************************************************

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nT, IsChi
real(kind=wp), intent(in) :: Zeta(nT), Eta(nT), P(nT,3), Q(nT,3), ChiI2
real(kind=wp), intent(out) :: T(nT), ZEInv(nT)
integer(kind=iwp) :: iT
real(kind=wp) :: PQ2

#include "macros.fh"
unused_var(Eta)
unused_var(IsChi)
unused_var(ChiI2)

#ifdef _DEBUGPRINT_
call RecPrt(' Zeta in TNAI1',' ',Zeta,nT,1)
call RecPrt(' Eta in TNAI1',' ',Eta,nT,1)
call RecPrt(' P in TNAI1',' ',P,nT,3)
call RecPrt(' Q in TNAI1',' ',Q,nT,3)
#endif
do iT=1,nT
  PQ2 = (P(iT,1)-Q(iT,1))**2+(P(iT,2)-Q(iT,2))**2+(P(iT,3)-Q(iT,3))**2
  T(iT) = Zeta(iT)*PQ2
end do
ZEInv(:) = One/Zeta

#ifdef _DEBUGPRINT_
call RecPrt('Tvalue',' ',T,nT,1)
#endif

return

end subroutine TNAI1
