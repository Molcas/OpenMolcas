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

subroutine t3sglh222(w,dima,dimb,dimc,s2,d2,ns)
! this routine adds following contribution to W
! for syma>symb;symc
!
! W(a,b;c) <- - S2 _i(b) . D2 _jk(a,c)
!
! w    - W matrix (I/O)
! dima - dimension of a index (I)
! dimb - dimension of b index (I)
! dimc - dimension of c index (I)
! s2   - S2 matrix (I)
! d2   - D2 matrix (I)
! ns   - signum of the contribution (+-1) (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dimc, ns
real(kind=wp), intent(inout) :: w(dima,dimb,dimc)
real(kind=wp), intent(in) :: s2(dimb), d2(dima,dimc)
integer(kind=iwp) :: a, b, c
real(kind=wp) :: s

if (ns == 1) then
  ! phase + 1

  do c=1,dimc
    do b=1,dimb
      s = s2(b)
      do a=1,dima
        w(a,b,c) = w(a,b,c)-d2(a,c)*s
      end do
    end do
  end do

else
  ! phase - 1

  do c=1,dimc
    do b=1,dimb
      s = s2(b)
      do a=1,dima
        w(a,b,c) = w(a,b,c)+d2(a,c)*s
      end do
    end do
  end do

end if

return

end subroutine t3sglh222
