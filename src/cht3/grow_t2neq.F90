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

implicit none
integer a, b, dima, dimb, nv, no, i, j
integer lasta, lastb
real*8 t2(nv,nv,no,no)
real*8 tmp(dima,dimb,no,no)

!mp write(6,*) 'grow_t2neq dima , dimb  ',dima,dimb
!mp write(6,*) 'grow_t2neq lasta, lastb ',lasta,lastb
!mp write(6,*) 'grow_t2neq no           ',no

!? if (lasta == lastb) then
!?   do j=1,no
!?     do i=1,no
!?       do a=1,dima
!?         do b=1,a
!?           t2(lasta+a,lastb+b,i,j) = 1.0d0*tmp(a,b,i,j)
!?cmp        if (a /= b) t2(lastb+b,lasta+a,j,i) = -1.0d0*tmp(a,b,j,i)
!?           if (a /= b) t2(lastb+b,lasta+a,j,i) = 1.0d0*tmp(a,b,i,j)
!?         end do
!?       end do
!?     end do
!?   end do

!? else

do j=1,no
  do i=1,no
    do b=1,dimb
      do a=1,dima
        t2(lasta+a,lastb+b,i,j) = 1.0d0*tmp(a,b,i,j)
        !mp t2(lastb+b,lasta+a,j,i) = -1.0d0*tmp(b,a,j,i)
        t2(lastb+b,lasta+a,j,i) = 1.0d0*tmp(a,b,i,j)
      end do
    end do
  end do
end do

!? end if

return

end subroutine grow_t2neq
