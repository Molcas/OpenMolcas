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

subroutine grow_t2_blocked(t2,tmp,dima,dimb,no,lasta,lastb,length1,length2,sym,switch)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dima, dimb, no, lasta, lastb, length1, length2
real(kind=wp) :: t2(length1,length2,no,no), tmp(dima,dimb,no,no)
logical(kind=iwp) :: sym, switch
integer(kind=iwp) :: a, b, i, j

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
      do a=1,dima
        if (.not. switch) then
          t2(lasta+a,lastb+b,i,j) = tmp(a,b,i,j)
        else
          !mp !t2(lasta+a,lastb+b,i,j) = tmp(a,b,j,i)
          t2(lasta+a,lastb+b,i,j) = tmp(a,b,i,j)
        end if
        !mpn t2(lastb+b,lasta+a,j,i) = tmp(a,b,i,j)

        if (sym) t2(lastb+b,lasta+a,j,i) = tmp(a,b,i,j)

      end do
    end do
  end do
end do

!? end if

return

end subroutine grow_t2_blocked
