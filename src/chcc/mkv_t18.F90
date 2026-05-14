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

subroutine MkV_T18(Va,V,dima,no)
! this routine does:
! Va(a,j,i,u) = - [2V(ai|ju)-V(aj|iu)]
! from V(o_a,P,Q,a)
!
! N.B. Kvajto odflaknute, ozaj ze hnus, treba sa zamysliet
! ci je nutva permutacia s a-ckom na konci - bod QK2.2 @@

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, no
real(kind=wp), intent(out) :: Va(dima,no,no,no)
real(kind=wp), intent(in) :: V(no,no,no,dima)
integer(kind=iwp) :: i, j, u

do u=1,no
  do i=1,no
    do j=1,no
      Va(:,j,i,u) = V(j,i,u,:)-Two*V(i,j,u,:)
    end do
  end do
end do

return

end subroutine MkV_T18
