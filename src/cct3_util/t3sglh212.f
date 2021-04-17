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
       subroutine t3sglh212 (w,dima,dimab,dimc,s1,d1,ns)
c
c     this routine add following contribution to W
c     for syma=symb;symc
c
c     W(ab,c) <-  + S1 _i(c) . D1 _jk(ab)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a (b) index (I)
c     dimab  - dimension of ab (ac,bc) index (I)
c     dimc   - dimension of c index (I)
c     s1     - S1 matrix (I)
c     d1     - D1 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimab,dimc,ns
       real*8 w(1:dimab,1:dimc)
       real*8 s1(1:dimc)
       real*8 d1(1:dimab)
c
c     help variables
c
       integer c,ab
       real*8 s
c
       if (ns.eq.1) then
c     phase + 1
c
       do 100 c=1,dimc
       s=s1(c)
       do 101 ab=1,dimab
       w(ab,c)=w(ab,c)+d1(ab)*s
 101    continue
 100    continue
c
       else
c     phase - 1
c
       do 200 c=1,dimc
       s=s1(c)
       do 201 ab=1,dimab
       w(ab,c)=w(ab,c)-d1(ab)*s
 201    continue
 200    continue
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(dima)
       end
