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

subroutine t3sglh323(w,dima,dimb,dimc,s1,d1,ns)
! this routine adds following contribution to W
! for syma;symb>symc
!
! W(a;b,c) <- + S1 _i(a) . D1 _jk(b,c)
!
! w    - W matrix (I/O)
! dima - dimension of a index (I)
! dimb - dimension of b index (I)
! dimc - dimension of c index (I)
! s1   - S1 matrix (I)
! d1   - D1 matrix (I)
! ns   - signum of the contribution (+-1) (I)

integer dima, dimb, dimc, ns
real*8 w(1:dima,1:dimb,1:dimc)
real*8 s1(1:dima)
real*8 d1(1:dimb,1:dimc)
! help variables
integer a, b, c
real*8 s

if (ns == 1) then
  ! phase +1

  do c=1,dimc
    do b=1,dimb
      s = d1(b,c)
      do a=1,dima
        w(a,b,c) = w(a,b,c)+s1(a)*s
      end do
    end do
  end do

else
  ! phase -1

  do c=1,dimc
    do b=1,dimb
      s = d1(b,c)
      do a=1,dima
        w(a,b,c) = w(a,b,c)-s1(a)*s
      end do
    end do
  end do

end if

return

end subroutine t3sglh323
