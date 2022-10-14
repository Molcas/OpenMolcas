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
       subroutine t3sglh312 (w,dima,dimb,dimbc,s1,d1,ns)
!
!     this routine add following contribution to W
!     for syma;symb=symc
!
!     W(a;bc)  <- + S1 _i(a) . D1 _jk(bc)
!
!     w      - W  matrix (I/O)
!     dima   - dimension of a index (I)
!     dimb   - dimension of b (c) index (I)
!     dimbc  - dimension of bc index (I)
!     s1     - S1 matrix (I)
!     d1     - D1 matrix (I)
!     ns     - signum of the contribution (+-1) (I)
!
       integer dima,dimb,dimbc,ns
       real*8 w(1:dima,1:dimbc)
       real*8 s1(1:dima)
       real*8 d1(1:dimbc)
!
!     help variables
!
       integer a,bc
       real*8 s
!
       if (ns.eq.1) then
!     phase +1
!
       do 100 bc=1,dimbc
       s=d1(bc)
       do 101 a=1,dima
       w(a,bc)=w(a,bc)+s1(a)*s
 101    continue
 100    continue
!
       else
!     phase - 1
!
       do 200 bc=1,dimbc
       s=d1(bc)
       do 201 a=1,dima
       w(a,bc)=w(a,bc)-s1(a)*s
 201    continue
 200    continue
!
       end if
!
       return
! Avoid unused argument warnings
      if (.false.) call Unused_integer(dimb)
       end
