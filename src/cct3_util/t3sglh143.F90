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
       subroutine t3sglh143 (w,dima,dimb,dimc,s3,d3,ns)
!
!     this routine add following contribution to W
!     for syma>symb>symc
!
!     W(a,b,c) <- + S3 _i(c) . D3 _jk(a,b)
!
!     w      - W  matrix (I/O)
!     dima   - dimension of a index (I)
!     dimb   - dimension of b index (I)
!     dimc   - dimension of c index (I)
!     s3     - S3 matrix (I)
!     d3     - D3 matrix (I)
!     ns     - signum of the contribution (+-1) (I)
!
       integer dima,dimb,dimc,ns
       real*8 w(1:dima,1:dimb,1:dimc)
       real*8 s3(1:dimc)
       real*8 d3(1:dima,1:dimb)
!
!     help variables
!
       integer a,b,c
       real*8 s
!
       if (ns.eq.1) then
!     phase +1
!
       do 120 c=1,dimc
       s=s3(c)
       do 121 b=1,dimb
       do 122 a=1,dima
       w(a,b,c)=w(a,b,c)+d3(a,b)*s
 122    continue
 121    continue
 120    continue
!
       else
!     phase = -1
!
       do 220 c=1,dimc
       s=s3(c)
       do 221 b=1,dimb
       do 222 a=1,dima
       w(a,b,c)=w(a,b,c)-d3(a,b)*s
 222    continue
 221    continue
 220    continue
!
       end if
!
       return
       end
