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

subroutine grow_vvoo_blocked(AA,BB,no,dima,dimb,lasta,lastb,length1,length2,sym)
! this routine does:
!
! grow A(1324)/(vvoo) from the blocked cholesky vectors B(1234)/(vovo)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: no, dima, dimb, lasta, lastb, length1, length2
real(kind=wp) :: AA(length1,length2,no,no), BB(dima,no,dimb,no)
logical(kind=iwp) :: sym
integer(kind=iwp) :: i1, i2, i3, i4

!mp write(u6,'(A,4(i3,x))') 'AA lasta, lastb, dima, dimb ',lasta,lastb,dima,dimb
!mp write(u6,'(A,2(i9,x))') 'chk_a ',lasta+dima,length1
!mp write(u6,'(A,2(i9,x))') 'chk_b ',lastb+dimb,length2

do i4=1,no
  do i3=1,no
    do i1=1,dima
      do i2=1,dimb
        AA(lasta+i1,lastb+i2,i3,i4) = BB(i1,i3,i2,i4)

        if (sym) AA(lastb+i2,lasta+i1,i4,i3) = BB(i1,i3,i2,i4)

      end do
    end do
  end do
end do

return

end subroutine grow_vvoo_blocked
