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

subroutine T2d_div(T2,OE,dima,dimb,adda,addb,no,nv)
! this routine does:
! T2(a',b',i,j) = T2(a',b',i,j/(e(i)+e(j)-e(a)-e(b))
! for aGrp=beGrp, where only a'>=b',i,j are valid,
! and completed also cases a'<b',i,j
!
! division of T2n amplitides by denominator

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, adda, addb, no, nv
real(kind=wp), intent(inout) :: T2(dima,dimb,no,no)
real(kind=wp), intent(in) :: OE(no+nv)
integer(kind=iwp) :: a, av, b, bv, i, j

av = no+adda
bv = no+addb

!1 division by denominators

do j=1,no
  do i=1,no
    do a=1,dima
      T2(a,1:a,i,j) = T2(a,1:a,i,j)/(OE(i)+OE(j)-OE(av+a)-OE(bv+1:bv+a))
    end do
  end do
end do

!2 completing upper triangle

do j=1,no
  do i=1,no
    do b=2,dima
      T2(1:b-1,b,i,j) = T2(b,1:b-1,j,i)
    end do
  end do
end do

return

end subroutine T2d_div
