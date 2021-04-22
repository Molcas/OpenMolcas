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
       subroutine t3aphlp7 (a2,a3,b,dimp,dimq,dimr,ns,szkey)
c
c     this routine realize following procedure
c     for symp.ne.symq.ne.symr
c
c     B(p,q,r)=B(p,q,r)+ns*(-A2(p,r,q)+A3(p,q,r))
c
c     a2     - A2 matrix (I)
c     a3     - A3 matrix (I)
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
       real*8 a2(1:dimp,1:dimr,1:dimq)
       real*8 a3(1:dimp,1:dimq,1:dimr)
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
       do 102 r=1,dimr
       do 1020 q=1,dimq
       do 1021 p=1,dimp
       b(p,q,r)=b(p,q,r)+a3(p,q,r)
 1021   continue
 1020   continue
 102    continue
c
       do 104 r=1,dimr
       do 1040 q=1,dimq
       do 1041 p=1,dimp
       b(p,q,r)=b(p,q,r)-a2(p,r,q)
 1041   continue
 1040   continue
 104    continue
c
       else
c     phase -1
c
       do 202 r=1,dimr
       do 2020 q=1,dimq
       do 2021 p=1,dimp
       b(p,q,r)=b(p,q,r)-a3(p,q,r)
 2021   continue
 2020   continue
 202    continue
c
       do 204 r=1,dimr
       do 2040 q=1,dimq
       do 2041 p=1,dimp
       b(p,q,r)=b(p,q,r)+a2(p,r,q)
 2041   continue
 2040   continue
 204    continue
c
       end if
c
       return
       end
