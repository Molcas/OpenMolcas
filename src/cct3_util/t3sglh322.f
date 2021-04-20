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
       subroutine t3sglh322 (w,dima,dimb,dimc,s2,d2,ns)
c
c     this routine add following contribution to W
c     for syma;symb>symc
c
c     W(a;b,c) <- - S2 _i(c) . D2 _jk(a,b)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a index (I)
c     dimb   - dimension of b index (I)
c     dimc   - dimension of c index (I)
c     s2     - S2 matrix (I)
c     d2     - D2 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimb,dimc,ns
       real*8 w(1:dima,1:dimb,1:dimc)
       real*8 s2(1:dimc)
       real*8 d2(1:dima,1:dimb)
c
c     help variables
c
       integer a,b,c
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       do 110 c=1,dimc
       s=s2(c)
       do 111 b=1,dimb
       do 112 a=1,dima
       w(a,b,c)=w(a,b,c)-d2(a,b)*s
 112    continue
 111    continue
 110    continue
c
       else
c     phase -1
c
       do 210 c=1,dimc
       s=s2(c)
       do 211 b=1,dimb
       do 212 a=1,dima
       w(a,b,c)=w(a,b,c)+d2(a,b)*s
 212    continue
 211    continue
 210    continue
c
       end if
c
       return
       end
