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
       subroutine t3sglh11 (w,dima,dimab,dimabc,s1,d1,ns)
c
c     this routine add following contribution to W
c     for syma=symb=symc
c
c     W(abc)  <-  + S1 _i(a) . D1 _jk(bc)
c     - S1 _i(b) . D1 _jk(ac)
c     + S1 _i(c) . D1 _jk(ab)
c
c     w      - W  matrix (I/O)
c     dima   - dimension of a (b,c) index (I)
c     dimab  - dimension of ab (ac,bc) index (I)
c     dimabc - dimension of abc index (I)
c     s1     - S1 matrix (I)
c     d1     - D1 matrix (I)
c     ns     - signum of the contribution (+-1) (I)
c
#include "t31.fh"
       integer dima,dimab,dimabc,ns
       real*8 w(1:dimabc)
       real*8 s1(1:dima)
       real*8 d1(1:dimab)
c
c     help variables
c
       integer a,b,c,ab0,ac0,bc0,abc
       real*8 s
c
       if (ns.eq.1) then
c     phase +1
c
       abc=0
       do 100 a=3,dima
       s=s1(a)
       do 101 b=2,a-1
       bc0=nshf(b)
       do 102 c=1,b-1
       abc=abc+1
       w(abc)=w(abc)+d1(bc0+c)*s
 102    continue
 101    continue
 100    continue
c
       abc=0
       do 110 a=3,dima
       ac0=nshf(a)
       do 111 b=2,a-1
       s=s1(b)
       do 112 c=1,b-1
       abc=abc+1
       w(abc)=w(abc)-d1(ac0+c)*s
 112    continue
 111    continue
 110    continue
c
       abc=0
       do 120 a=3,dima
       ab0=nshf(a)
       do 121 b=2,a-1
       s=d1(ab0+b)
       do 122 c=1,b-1
       abc=abc+1
       w(abc)=w(abc)+s1(c)*s
 122    continue
 121    continue
 120    continue
c
       else
c     phase -1
c
       abc=0
       do 200 a=3,dima
       s=s1(a)
       do 201 b=2,a-1
       bc0=nshf(b)
       do 202 c=1,b-1
       abc=abc+1
       w(abc)=w(abc)-d1(bc0+c)*s
 202    continue
 201    continue
 200    continue
c
       abc=0
       do 210 a=3,dima
       ac0=nshf(a)
       do 211 b=2,a-1
       s=s1(b)
       do 212 c=1,b-1
       abc=abc+1
       w(abc)=w(abc)+d1(ac0+c)*s
 212    continue
 211    continue
 210    continue
c
       abc=0
       do 220 a=3,dima
       ab0=nshf(a)
       do 221 b=2,a-1
       s=d1(ab0+b)
       do 222 c=1,b-1
       abc=abc+1
       w(abc)=w(abc)-s1(c)*s
 222    continue
 221    continue
 220    continue
c
       end if
c
       return
       end
