************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
       subroutine t3sglh122 (w,dima,dimab,dimc,s3,d3,ns)
c
c     this routine add following contribution to W
c     for syma=symb > symc
c
c     W(ab,c) <-  + S3 _i(c) . D3 _jk(ab)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a (b) index (I)
c     dimab   - dimension of ab index (I)
c     dimc   - dimension of c index (I)
c     s3     - S3 matrix (I)
c     d3     - D3 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimab,dimc,ns
       real*8 w(1:dimab,1:dimc)
       real*8 s3(1:dimc)
       real*8 d3(1:dimab)
c
c     help variables
c
       integer c,ab
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       do 120 c=1,dimc
       s=s3(c)
       do 121 ab=1,dimab
       w(ab,c)=w(ab,c)+d3(ab)*s
 121    continue
 120    continue
c
       else
c     phase -1
c
       do 220 c=1,dimc
       s=s3(c)
       do 221 ab=1,dimab
       w(ab,c)=w(ab,c)-d3(ab)*s
 221    continue
 220    continue
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(dima)
       end
