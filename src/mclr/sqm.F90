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

subroutine SqM(a,b,n)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: a(nTri_Elem(n))
real(kind=wp), intent(out) :: b(n,n)
integer(kind=iwp) :: i, ij

ij = 0
do i=1,n
  B(i,i:n) = A(ij+1:ij+n-i+1)
  B(i:n,i) = A(ij+1:ij+n-i+1)
  ij = ij+n-i+1
end do

return

end subroutine SqM
