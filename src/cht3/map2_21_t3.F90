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

subroutine Map2_21_t3(A,B,d1,d2)
! this routine does:
! map B(21) <- A(12)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: d1, d2
real(kind=wp), intent(in) :: A(d1,d2)
real(kind=wp), intent(out) :: B(d2,d1)
integer(kind=iwp) :: i1

do i1=1,d1
  B(:,i1) = A(i1,:)
end do

return

end subroutine Map2_21_t3
