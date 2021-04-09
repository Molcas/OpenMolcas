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
       subroutine t3sglh141 (w,dima,dimb,dimc,s1,d1,ns)
c
c     this routine add following contribution to W
c     for syma>symb>symc
c
c     W(a,b,c) <- + S1 _i(a) . D1 _jk(b,c)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a index (I)
c     dimb   - dimension of b index (I)
c     dimc   - dimension of c index (I)
c     s1     - S1 matrix (I)
c     d1     - D1 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimb,dimc,ns
       real*8 w(1:dima,1:dimb,1:dimc)
       real*8 s1(1:dima)
       real*8 d1(1:dimb,1:dimc)
c
c     help variables
c
       integer a,b,c
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       do 100 c=1,dimc
       do 101 b=1,dimb
       s=d1(b,c)
       do 102 a=1,dima
       w(a,b,c)=w(a,b,c)+s1(a)*s
 102    continue
 101    continue
 100    continue
c
       else
c     phase = -1
c
       do 200 c=1,dimc
       do 201 b=1,dimb
       s=d1(b,c)
       do 202 a=1,dima
       w(a,b,c)=w(a,b,c)-s1(a)*s
 202    continue
 201    continue
 200    continue
c
       end if
c
       return
       end
