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

subroutine print_ddNorm(nMOs,Wdd,det3)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nMOs
real(kind=wp), intent(in) :: Wdd(nMOs,nMOs), det3
integer(kind=iwp) :: i, j

do i=1,nMOs
  write(u6,'(20F8.3)') (Wdd(i,j)*det3,j=1,nMOs)
end do

return

end subroutine print_ddNorm
