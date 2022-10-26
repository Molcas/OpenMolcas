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

subroutine Map4_3412_t3(A,B,d1,d2,d3,d4)
! this routine does:
! map B(3412) <- A(1234)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: d1, d2, d3, d4
real(kind=wp), intent(in) :: A(d1,d2,d3,d4)
real(kind=wp), intent(out) :: B(d3,d4,d1,d2)
integer(kind=iwp) :: i2, i3, i4

do i2=1,d2
  do i4=1,d4
    do i3=1,d3
      B(i3,i4,:,i2) = A(:,i2,i3,i4)
    end do
  end do
end do

return

end subroutine Map4_3412_t3
