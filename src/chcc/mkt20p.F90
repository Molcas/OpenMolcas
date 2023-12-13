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

subroutine MkT20p(T2,V,oe,dima,adda,no)
! this routine does:
! T2(a'b',i,j)=(a'i|b'j)/Dija'b'
!  for aGrp=bGrp
! N.B.2 qvajt odflaknute, neurobene

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, adda, no
real(kind=wp), intent(out) :: T2(nTri_Elem(dima),no,no)
real(kind=wp), intent(in) :: V(dima,no,dima,no), oe(*)
integer(kind=iwp) :: a, ab, i, j

do j=1,no
  do i=1,no
    ab = 0
    do a=1,dima
      T2(ab+1:ab+a,i,j) = V(a,i,1:a,j)/(oe(i)+oe(j)-oe(adda+a)-oe(adda+1:adda+a))
      ab = ab+a
    end do
  end do
end do

return

end subroutine MkT20p
