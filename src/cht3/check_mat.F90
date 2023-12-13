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

subroutine check_mat(mat,dima,dimb)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: dima, dimb
real(kind=wp), intent(in) :: mat(dima,dimb)
integer(kind=iwp) :: i, j

do i=1,dima
  do j=1,dimb
    if (abs(mat(i,j)) > 1.0e5_wp) write(u6,*) 'i,j,mat(i,j) ',i,j,mat(i,j)
  end do
end do

return

end subroutine check_mat
