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
       subroutine t3aphlp5 (a1,a2,b,dimp,dimq,dimr,ns,szkey)
c
c     this routine realize following procedure
c     for symp.ne.symq.ne.symr
c
c     B(p,q,r)=B(p,q,r)+ns*(A1(q,r,p)-A2(p,r,q))
c
c     a1     - A1 matrix (I)
c     a2     - A2 matrix (I)
c     b      - B  matrix (I/O)
c     dimp   - dimension of p index (I)
c     dimq   - dimension of q index (I)
c     dimr   - dimension of r index (I)
c     ns     - singum of the permutation (+-1) (I)
c     szkey  - set zero key (I)
c     = 0 no vanishing
c     = 1 set B=0 at the beggining
c
       integer dimp,dimq,dimr,ns,szkey
       real*8 a1(1:dimq,1:dimr,1:dimp)
       real*8 a2(1:dimp,1:dimr,1:dimq)
       real*8 b(1:dimp,1:dimq,1:dimr)
c
c     help variables
c
       integer p,q,r
       integer nhelp
c
c
       if (szkey.eq.1) then
       nhelp=dimp*dimq*dimr
       call cct3_mv0zero (nhelp,nhelp,b)
       end if
c
       if (ns.eq.1) then
c     phase +1
c
       do 104 r=1,dimr
       do 1040 q=1,dimq
       do 1041 p=1,dimp
       b(p,q,r)=b(p,q,r)-a2(p,r,q)
 1041   continue
 1040   continue
 104    continue
c
       do 106 r=1,dimr
       do 1060 q=1,dimq
       do 1061 p=1,dimp
       b(p,q,r)=b(p,q,r)+a1(q,r,p)
 1061   continue
 1060   continue
 106    continue
c
       else
c     phase -1
c
       do 204 r=1,dimr
       do 2040 q=1,dimq
       do 2041 p=1,dimp
       b(p,q,r)=b(p,q,r)+a2(p,r,q)
 2041   continue
 2040   continue
 204    continue
c
       do 206 r=1,dimr
       do 2060 q=1,dimq
       do 2061 p=1,dimp
       b(p,q,r)=b(p,q,r)-a1(q,r,p)
 2061   continue
 2060   continue
 206    continue
c
       end if
c
       return
       end
