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

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nv, d2, dima, dimb, dimc, lasta, lastb, lastc
real(kind=wp), intent(inout) :: w3(nTri_Elem(nv),d2)
real(kind=wp), intent(in) :: AA(dima,dimb,dimc)
integer(kind=iwp) :: a, a_old, a_point, ab, b, b_old, b_point

if ((dima == 0) .or. (dimb == 0)) then
  write(u6,*) 'dima, dimb = ',dima,dimb
  write(u6,*) 'zle je'
  call abend()
end if

a_point = 0
b_point = 0
ab = 0
!mp write(u6,'(A,3(i5))') 'lasta, lastb, lastc = ',lasta,lastb,lastc
!mp write(u6,'(A,2(i5))') 'dima, dimb          = ',dima,dimb

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

        !write(u6,*) 'b, b_old = ',b,b_old
        if ((b /= b_old) .or. (b == max(1,lastb+1))) then
          !write(u6,*) 'wft'
          b_point = b_point+1
          b_old = b
        end if

        !mp if (lastc == 0) write(u6,'(A,5(i5))') 'ab, a, b, a_point, b_point = ',ab,a,b,a_point,b_point
        w3(ab,lastc+1:lastc+dimc) = AA(a_point,b_point,:)

      end if
    end if

  end do
end do

return

end subroutine grow_w3
