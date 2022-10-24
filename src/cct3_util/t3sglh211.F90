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

subroutine t3sglh211(w,dima,dimab,dimc,s1,d1,ns)
! this routine adds following contribution to W
! for syma=symb;symc
!
! W(ab,c) <-  + S1 _i(a) . D1 _jk(b,c)
! - S1 _i(b) . D1 _jk(a,c)
!
! w     - W matrix (I/O)
! dima  - dimension of a (b) index (I)
! dimab - dimension of ab  index (I)
! dimc  - dimension of c index (I)
! s1    - S1 matrix (I)
! d1    - D1 matrix (I)
! ns    - signum of the contribution (+-1) (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimab, dimc, ns
real(kind=wp), intent(inout) :: w(dimab,dimc)
real(kind=wp), intent(in) :: s1(dima), d1(dima,dimc)
integer(kind=iwp) :: a, ab, b, c
real(kind=wp) :: s

if (ns == 1) then
  ! phase + 1

  do c=1,dimc
    ab = 0
    do a=2,dima
      s = s1(a)
      do b=1,a-1
        ab = ab+1
        w(ab,c) = w(ab,c)+d1(b,c)*s
      end do
    end do
  end do

  do c=1,dimc
    ab = 0
    do a=2,dima
      s = d1(a,c)
      do b=1,a-1
        ab = ab+1
        w(ab,c) = w(ab,c)-s1(b)*s
      end do
    end do
  end do

else
  ! phase - 1

  do c=1,dimc
    ab = 0
    do a=2,dima
      s = s1(a)
      do b=1,a-1
        ab = ab+1
        w(ab,c) = w(ab,c)-d1(b,c)*s
      end do
    end do
  end do

  do c=1,dimc
    ab = 0
    do a=2,dima
      s = d1(a,c)
      do b=1,a-1
        ab = ab+1
        w(ab,c) = w(ab,c)+s1(b)*s
      end do
    end do
  end do

end if

return

end subroutine t3sglh211
