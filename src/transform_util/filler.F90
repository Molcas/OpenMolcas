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

subroutine Filler(N,M,A)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, M
real(kind=wp), intent(out) :: A(N,M)
integer(kind=iwp) :: i, j, k

k = 0
do i=1,N
  do j=1,M
    k = k+1
    A(i,j) = One*j+0.1_wp*i+0.001_wp*k
  end do
end do

return

end subroutine Filler
