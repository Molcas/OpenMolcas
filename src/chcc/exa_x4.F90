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

subroutine ExA_X4(A,Ap,no)
! this routine does:
! Ap(i,u,v) <- A(ii,u,v)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: no
real(kind=wp), intent(in) :: A(nTri_Elem(no),no**2)
real(kind=wp), intent(out) :: Ap(no,no**2)
integer(kind=iwp) :: i, ii

do i=1,no
  ii = nTri_Elem(i)
  Ap(i,:) = A(ii,:)
end do

return

end subroutine ExA_X4
