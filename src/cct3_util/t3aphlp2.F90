!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
       subroutine t3aphlp2 (a1,a2,a3,b,dimp,dimq,dimr,dimpq,ns,szkey)
!
!     this routine realize following procedure
!     for symp=symq.ne.symr
!
!     B(pq,r)=B(pq,r)+ns.(A1(q,r,p)-A2(p,r,q)+A3(pq,r))
!
!     a1     - A1 matrix (I)
!     a2     - A2 matrix (I)
!     a3     - A3 matrix (I)
!     b      - B  matrix (I/O)
!     dimp   - dimension of p index (I)
!     dimq   - dimension of q index (I)
!     dimr   - dimension of r index (I)
!     dimpq  - dimension of pq (I)
!     ns     - singum of the permutation (+-1) (I)
!     szkey  - set zero key (I)
!     = 0 no vanishing
!     = 1 set B=0 at the beggining
!
#include "t31.fh"
       integer dimp,dimq,dimr,dimpq,ns,szkey
       real*8 a1(1:dimq,1:dimr,1:dimp)
       real*8 a2(1:dimp,1:dimr,1:dimq)
       real*8 a3(1:dimpq,1:dimr)
       real*8 b(1:dimpq,1:dimr)
!
!     help variables
!
       integer p,q,r,pq,pq0
       integer nhelp
!
!
       if (szkey.eq.1) then
       nhelp=dimpq*dimr
       call cct3_mv0zero (nhelp,nhelp,b)
       end if
!
       if (ns.eq.1) then
!     phase +1
!
       do 102 r=1,dimr
       do 1020 pq=1,dimpq
       b(pq,r)=b(pq,r)+a3(pq,r)
 1020   continue
 102    continue
!
       do 104 r=1,dimr
       do 1040 p=2,dimp
       pq0=nshf(p)
       do 1041 q=1,p-1
       b(pq0+q,r)=b(pq0+q,r)-a2(p,r,q)
 1041   continue
 1040   continue
 104    continue
!
       do 106 r=1,dimr
       do 1060 p=2,dimp
       pq0=nshf(p)
       do 1061 q=1,p-1
       b(pq0+q,r)=b(pq0+q,r)+a1(q,r,p)
 1061   continue
 1060   continue
 106    continue
!
       else
!     phase -1
!
       do 202 r=1,dimr
       do 2020 pq=1,dimpq
       b(pq,r)=b(pq,r)-a3(pq,r)
 2020   continue
 202    continue
!
       do 204 r=1,dimr
       do 2040 p=2,dimp
       pq0=nshf(p)
       do 2041 q=1,p-1
       b(pq0+q,r)=b(pq0+q,r)+a2(p,r,q)
 2041   continue
 2040   continue
 204    continue
!
       do 206 r=1,dimr
       do 2060 p=2,dimp
       pq0=nshf(p)
       do 2061 q=1,p-1
       b(pq0+q,r)=b(pq0+q,r)-a1(q,r,p)
 2061   continue
 2060   continue
 206    continue
!
       end if
!
       return
       end
