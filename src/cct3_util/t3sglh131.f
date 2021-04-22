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
       subroutine t3sglh131 (w,dima,dimb,dimbc,s1,d1,ns)
c
c     this routine add following contribution to W
c     for syma=symb > symc
c
c     W(a,bc) <-  + S1 _i(a) . D1 _jk(bc)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a (b) index (I)
c     dimab   - dimension of ab index (I)
c     dimc   - dimension of c index (I)
c     s1     - S1 matrix (I)
c     d1     - D1 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
       integer dima,dimb,dimbc,ns
       real*8 w(1:dima,1:dimbc)
       real*8 s1(1:dima)
       real*8 d1(1:dimbc)
c
c     help variables
c
       integer a,bc
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       do 100 bc=1,dimbc
       s=d1(bc)
       do 101 a=1,dima
       w(a,bc)=w(a,bc)+s1(a)*s
 101    continue
 100    continue
c
       else
c     phase -1
c
       do 200 bc=1,dimbc
       s=d1(bc)
       do 201 a=1,dima
       w(a,bc)=w(a,bc)-s1(a)*s
 201    continue
 200    continue
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(dimb)
       end
