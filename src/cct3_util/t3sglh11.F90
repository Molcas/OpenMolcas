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

#include "t31.fh"
integer dima, dimab, dimabc, ns
real*8 w(1:dimabc)
real*8 s1(1:dima)
real*8 d1(1:dimab)
! help variables
integer a, b, c, ab0, ac0, bc0, abc
real*8 s

if (ns == 1) then
  ! phase +1

  abc = 0
  do a=3,dima
    s = s1(a)
    do b=2,a-1
      bc0 = nshf(b)
      do c=1,b-1
        abc = abc+1
        w(abc) = w(abc)+d1(bc0+c)*s
      end do
    end do
  end do

  abc = 0
  do a=3,dima
    ac0 = nshf(a)
    do b=2,a-1
      s = s1(b)
      do c=1,b-1
        abc = abc+1
        w(abc) = w(abc)-d1(ac0+c)*s
      end do
    end do
  end do

  abc = 0
  do a=3,dima
    ab0 = nshf(a)
    do b=2,a-1
      s = d1(ab0+b)
      do c=1,b-1
        abc = abc+1
        w(abc) = w(abc)+s1(c)*s
      end do
    end do
  end do

else
  ! phase -1

  abc = 0
  do a=3,dima
    s = s1(a)
    do b=2,a-1
      bc0 = nshf(b)
      do c=1,b-1
        abc = abc+1
        w(abc) = w(abc)-d1(bc0+c)*s
      end do
    end do
  end do

  abc = 0
  do a=3,dima
    ac0 = nshf(a)
    do b=2,a-1
      s = s1(b)
      do c=1,b-1
        abc = abc+1
        w(abc) = w(abc)+d1(ac0+c)*s
      end do
    end do
  end do

  abc = 0
  do a=3,dima
    ab0 = nshf(a)
    do b=2,a-1
      s = d1(ab0+b)
      do c=1,b-1
        abc = abc+1
        w(abc) = w(abc)-s1(c)*s
      end do
    end do
  end do

end if

return

end subroutine t3sglh11
