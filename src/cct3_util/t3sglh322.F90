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

subroutine t3sglh322(w,dima,dimb,dimc,s2,d2,ns)
! this routine adds following contribution to W
! for syma;symb>symc
!
! W(a;b,c) <- - S2 _i(c) . D2 _jk(a,b)
!
! w    - W matrix (I/O)
! dima - dimension of a index (I)
! dimb - dimension of b index (I)
! dimc - dimension of c index (I)
! s2   - S2 matrix (I)
! d2   - D2 matrix (I)
! ns   - signum of the contribution (+-1) (I)

integer dima, dimb, dimc, ns
real*8 w(1:dima,1:dimb,1:dimc)
real*8 s2(1:dimc)
real*8 d2(1:dima,1:dimb)
! help variables
integer a, b, c
real*8 s

if (ns == 1) then
  ! phase +1

  do c=1,dimc
    s = s2(c)
    do b=1,dimb
      do a=1,dima
        w(a,b,c) = w(a,b,c)-d2(a,b)*s
      end do
    end do
  end do

else
  ! phase -1

  do c=1,dimc
    s = s2(c)
    do b=1,dimb
      do a=1,dima
        w(a,b,c) = w(a,b,c)+d2(a,b)*s
      end do
    end do
  end do

end if

return

end subroutine t3sglh322
