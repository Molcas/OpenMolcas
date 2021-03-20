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
! Copyright (C) Roland Lindh                                           *
!               Valera Veryazov                                        *
!***********************************************************************

subroutine save_ddNorm(ddNorm,ii,jj,Wdd,nMOs)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ii, jj, nMOs
real(kind=wp), intent(in) :: ddNorm
real(kind=wp), intent(inout) :: Wdd(nMOs,nMOs)

Wdd(ii,jj) = Wdd(ii,jj)+ddNorm

return

end subroutine save_ddNorm
