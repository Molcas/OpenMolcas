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

subroutine t3sglh141(w,dima,dimb,dimc,s1,d1,ns)
! this routine adds following contribution to W
! for syma>symb>symc
!
! W(a,b,c) <- + S1 _i(a) . D1 _jk(b,c)
!
! w    - W matrix (I/O)
! dima - dimension of a index (I)
! dimb - dimension of b index (I)
! dimc - dimension of c index (I)
! s1   - S1 matrix (I)
! d1   - D1 matrix (I)
! ns   - signum of the contribution (+-1) (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dimc, ns
real(kind=wp), intent(inout) :: w(dima,dimb,dimc)
real(kind=wp), intent(in) :: s1(dima), d1(dimb,dimc)
integer(kind=iwp) :: a

if (ns == 1) then
  ! phase +1

  do a=1,dima
    w(a,:,:) = w(a,:,:)+s1(a)*d1
  end do

else
  ! phase = -1

  do a=1,dima
    w(a,:,:) = w(a,:,:)-s1(a)*d1
  end do

end if

return

end subroutine t3sglh141
