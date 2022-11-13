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

subroutine t3sglh132(w,dima,dimb,dimbc,s2,d2,ns)
! this routine adds following contribution to W
! for syma=symb > symc
!
! W(a,bc) <-  - S2 _i(b) . D2 _jk(a,c)
!             + S2 _i(c) . D2 _jk(a,b)
!
! w     - W matrix (I/O)
! dima  - dimension of a (b) index (I)
! dimab - dimension of ab index (I)
! dimc  - dimension of c index (I)
! s2    - S2 matrix (I)
! d2    - D2 matrix (I)
! ns    - signum of the contribution (+-1) (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dimbc, ns
real(kind=wp), intent(inout) :: w(dima,dimbc)
real(kind=wp), intent(in) :: s2(dimb), d2(dima,dimb)
integer(kind=iwp) :: b, bc, c

if (ns == 1) then
  ! phase +1

  bc = 0
  do b=2,dimb
    w(:,bc+1:bc+b-1) = w(:,bc+1:bc+b-1)-d2(:,1:b-1)*s2(b)
    bc = bc+b-1
  end do

  bc = 0
  do b=2,dimb
    do c=1,b-1
      bc = bc+1
      w(:,bc) = w(:,bc)+d2(:,b)*s2(c)
    end do
  end do

else
  ! phase -1

  bc = 0
  do b=2,dimb
    w(:,bc+1:bc+b-1) = w(:,bc+1:bc+b-1)+d2(:,1:b-1)*s2(b)
    bc = bc+b-1
  end do

  bc = 0
  do b=2,dimb
    do c=1,b-1
      bc = bc+1
      w(:,bc) = w(:,bc)-d2(:,b)*s2(c)
    end do
  end do

end if

return

end subroutine t3sglh132
