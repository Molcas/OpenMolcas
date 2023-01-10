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

subroutine grow_l1(l1,tmp,dima,nc,no,nv,last)
! this routine does:
!
! grow Cholesky vectors L1(m,i,a) by the segment in tmp

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, nc, no, nv, last
real(kind=wp), intent(inout) :: l1(nc,no,nv)
real(kind=wp), intent(in) :: tmp(nc,no,dima)

!mp write(u6,*) 'grow_l1i ',dima

l1(:,:,last+1:last+dima) = tmp

return

end subroutine grow_l1
