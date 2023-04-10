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

implicit none
integer dima, dimb, adda, addb, no, nv
real*8 T2(1:dima,1:dimb,1:no,1:no)
real*8 OE(1:no+nv)
! help variables
integer i, j, a, av, b, bv
real*8 eijb

av = no+adda
bv = no+addb

do j=1,no
  do i=1,no
    do b=1,dimb
      eijb = OE(i)+OE(j)-OE(bv+b)
      do a=1,dima
        T2(a,b,i,j) = T2(a,b,i,j)/(eijb-OE(av+a))
      end do
    end do
  end do
end do

return

end subroutine T2od_div
