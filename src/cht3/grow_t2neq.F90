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

subroutine grow_t2neq(t2,tmp,dima,dimb,nv,no,lasta,lastb)
! this routine does:
!
! grow amplitude file t2(a,b,i,j) by the segment in tmp
! for case sa != sb

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, nv, no, lasta, lastb
real(kind=wp), intent(inout) :: t2(nv,nv,no,no)
real(kind=wp), intent(in) :: tmp(dima,dimb,no,no)
integer(kind=iwp) :: b, i, j

!mp write(u6,*) 'grow_t2neq dima , dimb  ',dima,dimb
!mp write(u6,*) 'grow_t2neq lasta, lastb ',lasta,lastb
!mp write(u6,*) 'grow_t2neq no           ',no

!? if (lasta == lastb) then
!?   do j=1,no
!?     do i=1,no
!?       do a=1,dima
!?         do b=1,a
!?           t2(lasta+a,lastb+b,i,j) = tmp(a,b,i,j)
!?cmp        if (a /= b) t2(lastb+b,lasta+a,j,i) = -tmp(a,b,j,i)
!?           if (a /= b) t2(lastb+b,lasta+a,j,i) = tmp(a,b,i,j)
!?         end do
!?       end do
!?     end do
!?   end do

!? else

do j=1,no
  do i=1,no
    do b=1,dimb
      t2(lasta+1:lasta+dima,lastb+b,i,j) = tmp(:,b,i,j)
      !mp t2(lastb+b,lasta+1:lasta+dima,j,i) = -tmp(b,1:dima,j,i)
      t2(lastb+b,lasta+1:lasta+dima,j,i) = tmp(:,b,i,j)
    end do
  end do
end do

!? end if

return

end subroutine grow_t2neq
