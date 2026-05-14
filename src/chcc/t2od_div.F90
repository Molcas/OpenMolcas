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

subroutine T2od_div(T2,OE,dima,dimb,adda,addb,no,nv)
! this routine does:
! T2(a',b',i,j) = T2(a',b',i,j/(e(i)+e(j)-e(a)-e(b))
!
! division of T1n amplitides by denominator

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, adda, addb, no, nv
real(kind=wp), intent(inout) :: T2(dima,dimb,no,no)
real(kind=wp), intent(in) :: OE(no+nv)
integer(kind=iwp) :: av, b, bv, i, j

av = no+adda
bv = no+addb

do j=1,no
  do i=1,no
    do b=1,dimb
      T2(:,b,i,j) = T2(:,b,i,j)/(OE(i)+OE(j)-OE(bv+b)-OE(av+1:av+dima))
    end do
  end do
end do

return

end subroutine T2od_div
