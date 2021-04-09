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
       subroutine exth5 (a,b,dimp,dimq,dimqr,q)
c
c     this routine extract A(p,qr) -> B_q(p,r)
c
c     a     - matrxi a (Input)
c     b     - matrix b (Output)
c     dimp  - dimension of p (Input)
c     dimq  - dimension of q (and also r) (Input)
c     dimqr - dimension of qr (Input)
c     q     - value of index q (Input)
c
#include "t31.fh"
       integer dimp,dimq,dimqr,q
       real*8 a(1:dimp,1:dimqr)
       real*8 b(1:dimp,dimq)
c
c     help variables
c
       integer r,p,rq,qr0,qr
c
       if (q.eq.0) then
       return
       end if
c
c     r>q part
       if (q.gt.1) then
       qr0=nshf(q)
       do 20 r=1,q-1
       qr=qr0+r
       do 21 p=1,dimp
       b(p,r)=a(p,qr)
 21     continue
 20     continue
       end if
c
c     r=q part
       do 40 p=1,dimp
       b(p,q)=0.0d0
 40     continue
c
c     r<p part
       if (q.lt.dimq) then
       do 60 r=q+1,dimq
       rq=nshf(r)+q
       do 61 p=1,dimp
       b(p,r)=-a(p,rq)
 61     continue
 60     continue
       end if
c
       return
       end
