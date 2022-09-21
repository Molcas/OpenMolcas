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

subroutine grow_w3(w3,AA,nv,d2,dima,dimb,dimc,lasta,lastb,lastc)
! this routine does:
!
! add the block contribution AA(a',b',c') to w3(a>=b,c)

implicit none
integer a, b, c, dima, dimb, dimc, lasta, lastb, lastc, ab, nv
integer a_point, b_point
integer d2
real*8 w3(nv*(nv+1)/2,d2)
real*8 AA(dima,dimb,dimc)
integer a_old, b_old

if ((dima == 0) .or. (dimb == 0)) then
  write(6,*) 'dima, dimb = ',dima,dimb
  write(6,*) 'zle je'
  call abend()
end if

a_point = 0
b_point = 0
ab = 0
!mp write(6,'(A,3(i5))') 'lasta, lastb, lastc = ',lasta,lastb,lastc
!mp write(6,'(A,2(i5))') 'dima, dimb          = ',dima,dimb

a_old = 0
b_old = 0

do a=1,nv
  b_point = 0
  do b=1,a
    ab = ab+1
    if ((a >= lasta+1) .and. (a <= lasta+dima)) then

      if (a /= a_old) then
        a_point = a_point+1
        a_old = a
      end if

      if ((b >= max(1,lastb+1)) .and. (b <= min(a,lastb+dimb))) then

        !write(6,*) 'b, b_old = ',b,b_old
        if ((b /= b_old) .or. (b == max(1,lastb+1))) then
          !write(6,*) 'wft'
          b_point = b_point+1
          b_old = b
        end if

        !mp if (lastc == 0) write(6,'(A,5(i5))') 'ab, a, b, a_point, b_point = ',ab,a,b,a_point,b_point
        do c=1,dimc
          w3(ab,lastc+c) = AA(a_point,b_point,c)
        end do

      end if
    end if

  end do
end do

return

end subroutine grow_w3
