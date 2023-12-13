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

subroutine grow_t2anti_blocked1(t2,tmp,dima,dimb,no,lasta,lastb,length1,length2)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, no, lasta, lastb, length1, length2
real(kind=wp), intent(inout) :: t2(length1,length2,nTri_Elem(no-1))
real(kind=wp), intent(in) :: tmp(dima,dimb,no,no)
integer(kind=iwp) :: i, ij, j

!mp write(u6,*) 'lasta+dima, length1 ',lasta+dima,length1
!mp write(u6,*) 'lastb+dimb, length2 ',lastb+dimb,length2

!mp write(u6,'(A,2(i10,x),i3)') 'length1, length2, no ',length1,length2,no
!mp write(u6,'(A,4(i3,x))') 'dima, dimb, lasta, lastb',dima,dimb,lasta,lastb

ij = 0
do i=2,no
  do j=1,i-1
    ij = ij+1
    t2(lasta+1:lasta+dima,lastb+1:lastb+dimb,ij) = tmp(:,:,i,j)-tmp(:,:,j,i)
  end do
end do

return

end subroutine grow_t2anti_blocked1
