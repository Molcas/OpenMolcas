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

subroutine grow_t2_blocked(t2,tmp,dima,dimb,no,lasta,lastb,length1,length2,sym)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, no, lasta, lastb, length1, length2
real(kind=wp), intent(inout) :: t2(length1,length2,no,no)
real(kind=wp), intent(in) :: tmp(dima,dimb,no,no)
logical(kind=iwp) :: sym
integer(kind=iwp) :: b, i, j

!mp write(u6,*) 'grow_t2neq dima , dimb  ',dima,dimb
!mp write(u6,*) 'grow_t2neq lasta, lastb ',lasta,lastb
!mp write(u6,*) 'grow_t2neq no           ',no

!mp if (lasta == lastb) then
!? if (grpa == grpb) then
!?   do j=1,no
!?     do i=1,no
!?       do a=1,dima
!?         do b=1,a
!?           t2(lasta+a,lastb+b,i,j) = tmp(a,b,i,j)
!?           if (a /= b) t2(lastb+b,lasta+a,j,i) = tmp(a,b,i,j)
!?         end do
!?       end do
!?     end do
!?   end do

!? else

do j=1,no
  do i=1,no
    do b=1,dimb
      !if (switch) then
      !  !mp !t2(lasta+1:lasta+dima,lastb+b,i,j) = tmp(:,b,j,i)
      !  t2(lasta+1:lasta+dima,lastb+b,i,j) = tmp(:,b,i,j)
      !else
      t2(lasta+1:lasta+dima,lastb+b,i,j) = tmp(:,b,i,j)
      !end if
      !mpn t2(lastb+b,lasta+a,j,i) = tmp(a,b,i,j)

      if (sym) t2(lastb+b,lasta+1:lasta+dima,j,i) = tmp(:,b,i,j)
    end do
  end do
end do

!? end if

return

end subroutine grow_t2_blocked
