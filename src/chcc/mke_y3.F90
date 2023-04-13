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

subroutine MkE_Y3(Va,V,dima,dimb,no)
! this routine does:
! Va(a',i,b',j) = 2 V(a',j,b',i) - V(a',i,b',j)
!
! N.B. Kvajt odflaknute

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dima, dimb, no
real(kind=wp) :: Va(dima,no,dimb,no), V(dima,no,dimb,no)
integer(kind=iwp) :: a, b, i, j

do j=1,no
  do b=1,dimb
    do i=1,no
      do a=1,dima
        Va(a,i,b,j) = Two*V(a,j,b,i)-V(a,i,b,j)
      end do
    end do
  end do
end do

return

end subroutine MkE_Y3
