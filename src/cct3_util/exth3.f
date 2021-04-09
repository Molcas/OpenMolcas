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
       subroutine exth3 (a,b,dimp,dimq,dimr,q,nfact)
c
c     this routine extract A(p,q,r) -> B_q(p,r)
c
c     a     - matrxi a (Input)
c     b     - matrix b (Output)
c     dimp  - dimension of p (Input)
c     dimq  - dimension of q (Input)
c     dimr  - dimension of r (Input)
c     q     - value of index q (Input)
c     nfact - sign (+-1,0) (Input)
c
       integer dimp,dimq,dimr,q,nfact
       real*8 a(1:dimp,1:dimq,1:dimr)
       real*8 b(1:dimp,1:dimr)
c
c     help variables
c
       integer p,r
c
       if (nfact.eq.1) then
       do 10 r=1,dimr
       do 11 p=1,dimp
       b(p,r)=a(p,q,r)
 11     continue
 10     continue
       else if (nfact.eq.-1) then
       do 20 r=1,dimr
       do 21 p=1,dimp
       b(p,r)=-a(p,q,r)
 21     continue
 20     continue
       else if (nfact.eq.0) then
       do 30 r=1,dimr
       do 31 p=1,dimp
       b(p,r)=0.0d0
 31     continue
 30     continue
       end if

c
       return
       end
