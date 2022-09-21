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

subroutine grow_t2_fblocked1(t2,tmp,dima,dimb,no,lasta,lastb,length1,length2)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dima, dimb, no, lasta, lastb, length1, length2
real(kind=wp) :: t2(length1,length2,no,no), tmp(dima,dimb,no,no)
integer(kind=iwp) :: a, b, i, j

!mp write(u6,*) 'lasta+dima, length1 ',lasta+dima,length1
!mp write(u6,*) 'lastb+dimb, length2 ',lastb+dimb,length2

!mp write(u6,'(A,2(i10,x),i3)') 'length1, length2, no ',length1,length2,no
!mp write(u6,'(A,4(i3,x))') 'dima, dimb, lasta, lastb',dima,dimb,lasta,lastb

!mp ij = 0
!mp do i=2,no
!mp   do j=1,i-1
!mp     ij = ij+1
do i=1,no
  do j=1,no
    do b=1,dimb
      do a=1,dima

        !mp t2(lasta+a,lastb+b,ij) = tmp(a,b,i,j)-tmp(a,b,j,i)
        t2(lasta+a,lastb+b,i,j) = tmp(a,b,i,j)

      end do
    end do
  end do
end do

return

end subroutine grow_t2_fblocked1
