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

subroutine exMap3_231(A,B,d1,d2)
! this routine does:
!
! A (a,bc) -> B(b,c,a)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: d1, d2
real(kind=wp), intent(in) :: A(d1,nTri_Elem(d2))
real(kind=wp), intent(out) :: B(d2,d2,d1)
integer(kind=iwp) :: i2, i23, i3

i23 = 0
do i2=1,d2
  do i3=1,i2
    i23 = i23+1
    B(i2,i3,:) = A(:,i23)
    B(i3,i2,:) = A(:,i23)
  end do
end do

return

end subroutine exMap3_231
