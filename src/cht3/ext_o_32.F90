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

subroutine ext_o_32(A,B,nc,no,dima,occ_ind)
! this routine does:
!
! extract B (m,a')_i <- A (m,i,a')

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nc, no, dima, occ_ind
real(kind=wp), intent(in) :: A(nc,no,dima)
real(kind=wp), intent(out) :: B(nc,dima)
integer(kind=iwp) :: i2

do i2=1,dima
  B(:,i2) = A(:,occ_ind,i2)
end do

return

end subroutine ext_o_32
