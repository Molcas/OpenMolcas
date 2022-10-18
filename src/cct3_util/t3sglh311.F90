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

subroutine t3sglh311(w,dima,dimb,dimbc,s1,d1,ns)
! this routine adds following contribution to W
! for syma;symb=symc
!
! W(a;bc)  <- + S1 _i(b) . D1 _jk(a,c)
! - S1 _i(c) . D2 _jk(a,b)
!
! w     - W matrix (I/O)
! dima  - dimension of a index (I)
! dimb  - dimension of b (c) index (I)
! dimbc - dimension of bc index (I)
! s1    - S1 matrix (I)
! s2    - S2 matrix (I)
! d1    - D1 matrix (I)
! d2    - D2 matrix (I)
! ns    - signum of the contribution (+-1) (I)

integer dima, dimb, dimbc, ns
real*8 w(1:dima,1:dimbc)
real*8 s1(1:dimb)
real*8 d1(1:dima,1:dimb)
! help variables
integer a, b, c, bc
real*8 s

if (ns == 1) then
  ! phase +1

  bc = 0
  do b=2,dimb
    s = s1(b)
    do c=1,b-1
      bc = bc+1
      do a=1,dima
        w(a,bc) = w(a,bc)+d1(a,c)*s
      end do
    end do
  end do

  bc = 0
  do b=2,dimb
    do c=1,b-1
      bc = bc+1
      s = s1(c)
      do a=1,dima
        w(a,bc) = w(a,bc)-d1(a,b)*s
      end do
    end do
  end do

else
  ! phase - 1

  bc = 0
  do b=2,dimb
    s = s1(b)
    do c=1,b-1
      bc = bc+1
      do a=1,dima
        w(a,bc) = w(a,bc)-d1(a,c)*s
      end do
    end do
  end do

  bc = 0
  do b=2,dimb
    do c=1,b-1
      bc = bc+1
      s = s1(c)
      do a=1,dima
        w(a,bc) = w(a,bc)+d1(a,b)*s
      end do
    end do
  end do

end if

return

end subroutine t3sglh311
