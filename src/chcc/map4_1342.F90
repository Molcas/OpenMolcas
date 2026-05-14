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

subroutine Map4_1342(A,B,d1,d2,d3,d4)
! this routine does:
! map B(1423) <- A(1234)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: d1, d2, d3, d4
real(kind=wp), intent(in) :: A(d1,d2,d3,d4)
real(kind=wp), intent(out) :: B(d1,d4,d2,d3)
integer(kind=iwp) :: i2, i3

do i3=1,d3
  do i2=1,d2
    b(:,:,i2,i3) = a(:,i2,i3,:)
  end do
end do

return

end subroutine Map4_1342
