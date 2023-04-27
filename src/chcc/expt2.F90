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

subroutine ExpT2(T2p,T2u,dima,dimab,no)
! this routine does:
! Make T2u(a',b',i,j) from T2p(a'b',i,j)
!
! N.B. Kvajt odflaknute, vypocet ab sa da dat zefektivnit

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimab, no
real(kind=wp), intent(in) :: T2p(dimab,no,no)
real(kind=wp), intent(out) :: T2u(dima,dima,no,no)
integer(kind=iwp) :: a, ab, b, ba0, i, j

do j=1,no
  do i=1,no
    do b=1,dima
      ba0 = nTri_Elem(b-1)
      T2u(1:b,b,i,j) = T2p(ba0+1:ba0+b,j,i)
      do a=1+b,dima
        ab = nTri_Elem(a-1)+b
        T2u(a,b,i,j) = T2p(ab,i,j)
      end do
    end do
  end do
end do

return

end subroutine ExpT2
