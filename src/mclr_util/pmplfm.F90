!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine PMPLFM(AP,B,NDIM)
! Add lower half of a full matrix to a matrix packed
! in lower triangular form (packed matrix stored columnwise)

use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: AP(*)
real(kind=wp), intent(in) :: B(*)
integer(kind=iwp), intent(in) :: NDIM
integer(kind=iwp) :: IBSF, IBSP, ICOL, NELMNT

IBSP = 1
IBSF = 1
do ICOL=1,NDIM
  NELMNT = NDIM-ICOL+1
  AP(IBSP:IBSP+NELMNT-1) = AP(IBSP:IBSP+NELMNT-1)+B(IBSF:IBSF+NELMNT-1)
  IBSP = IBSP+NELMNT
  IBSF = IBSF+NDIM
end do

return

end subroutine PMPLFM
