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

subroutine expand4_12(AA,BB,d1,d2,d3)
! this routine does:
!
! A(ab,i,j) -> A(a,b,i,j)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: d1, d2, d3
real(kind=wp), intent(in) :: AA(nTri_Elem(d1),d2,d3)
real(kind=wp), intent(out) :: BB(d1,d1,d2,d3)
integer(kind=iwp) :: a, ab, b, i

ab = 0
do a=1,d1
  do b=1,a-1
    ab = ab+1
    BB(a,b,:,:) = AA(ab,:,:)
    do i=1,d2
      BB(b,a,:,i) = AA(ab,i,:)
    end do
  end do
  ab = ab+1
  BB(a,a,:,:) = AA(ab,:,:)
end do

return

end subroutine expand4_12
