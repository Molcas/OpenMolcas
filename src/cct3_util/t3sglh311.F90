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
       subroutine t3sglh311 (w,dima,dimb,dimbc,s1,d1,ns)
!
!     this routine add following contribution to W
!     for syma;symb=symc
!
!     W(a;bc)  <- + S1 _i(b) . D1 _jk(a,c)
!     - S1 _i(c) . D2 _jk(a,b)
!
!     w      - W  matrix (I/O)
!     dima   - dimension of a index (I)
!     dimb   - dimension of b (c) index (I)
!     dimbc  - dimension of bc index (I)
!     s1     - S1 matrix (I)
!     s2     - S2 matrix (I)
!     d1     - D1 matrix (I)
!     d2     - D2 matrix (I)
!     ns     - signum of the contribution (+-1) (I)
!
       integer dima,dimb,dimbc,ns
       real*8 w(1:dima,1:dimbc)
       real*8 s1(1:dimb)
       real*8 d1(1:dima,1:dimb)
!
!     help variables
!
       integer a,b,c,bc
       real*8 s
!
       if (ns.eq.1) then
!     phase +1
!
       bc=0
       do 100 b=2,dimb
       s=s1(b)
       do 101 c=1,b-1
       bc=bc+1
       do 102 a=1,dima
       w(a,bc)=w(a,bc)+d1(a,c)*s
 102    continue
 101    continue
 100    continue
!
       bc=0
       do 110 b=2,dimb
       do 111 c=1,b-1
       bc=bc+1
       s=s1(c)
       do 112 a=1,dima
       w(a,bc)=w(a,bc)-d1(a,b)*s
 112    continue
 111    continue
 110    continue
!
       else
!     phase - 1
!
       bc=0
       do 200 b=2,dimb
       s=s1(b)
       do 201 c=1,b-1
       bc=bc+1
       do 202 a=1,dima
       w(a,bc)=w(a,bc)-d1(a,c)*s
 202    continue
 201    continue
 200    continue
!
       bc=0
       do 210 b=2,dimb
       do 211 c=1,b-1
       bc=bc+1
       s=s1(c)
       do 212 a=1,dima
       w(a,bc)=w(a,bc)+d1(a,b)*s
 212    continue
 211    continue
 210    continue
!
       end if
!
       return
       end
