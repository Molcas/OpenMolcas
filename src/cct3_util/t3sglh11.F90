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

subroutine t3sglh11(w,dima,dimab,dimabc,s1,d1,ns)
! this routine adds following contribution to W
! for syma=symb=symc
!
! W(abc)  <-  + S1 _i(a) . D1 _jk(bc)
! - S1 _i(b) . D1 _jk(ac)
! + S1 _i(c) . D1 _jk(ab)
!
! w      - W matrix (I/O)
! dima   - dimension of a (b,c) index (I)
! dimab  - dimension of ab (ac,bc) index (I)
! dimabc - dimension of abc index (I)
! s1     - S1 matrix (I)
! d1     - D1 matrix (I)
! ns     - signum of the contribution (+-1) (I)

use CCT3_global, only: nshf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimab, dimabc, ns
real(kind=wp), intent(inout) :: w(dimabc)
real(kind=wp), intent(in) :: s1(dima), d1(dimab)
integer(kind=iwp) :: a, ab0, abc, ac0, b, bc0

if (ns == 1) then
  ! phase +1

  abc = 0
  do a=3,dima
    do b=2,a-1
      bc0 = nshf(b)
      w(abc+1:abc+b-1) = w(abc+1:abc+b-1)+d1(bc0+1:bc0+b-1)*s1(a)
      abc = abc+b-1
    end do
  end do

  abc = 0
  do a=3,dima
    ac0 = nshf(a)
    do b=2,a-1
      w(abc+1:abc+b-1) = w(abc+1:abc+b-1)-d1(ac0+1:ac0+b-1)*s1(b)
      abc = abc+b-1
    end do
  end do

  abc = 0
  do a=3,dima
    ab0 = nshf(a)
    do b=2,a-1
      w(abc+1:abc+b-1) = w(abc+1:abc+b-1)+s1(1:b-1)*d1(ab0+b)
      abc = abc+b-1
    end do
  end do

else
  ! phase -1

  abc = 0
  do a=3,dima
    do b=2,a-1
      bc0 = nshf(b)
      w(abc+1:abc+b-1) = w(abc+1:abc+b-1)-d1(bc0+1:bc0+b-1)*s1(a)
      abc = abc+b-1
    end do
  end do

  abc = 0
  do a=3,dima
    ac0 = nshf(a)
    do b=2,a-1
      w(abc+1:abc+b-1) = w(abc+1:abc+b-1)+d1(ac0+1:ac0+b-1)*s1(b)
      abc = abc+b-1
    end do
  end do

  abc = 0
  do a=3,dima
    ab0 = nshf(a)
    do b=2,a-1
      w(abc+1:abc+b-1) = w(abc+1:abc+b-1)-s1(1:b-1)*d1(ab0+b)
      abc = abc+b-1
    end do
  end do

end if

return

end subroutine t3sglh11
