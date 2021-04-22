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
       subroutine t3dhlp1 (w,v,dimp,dimq,dimr,denijk,ec,
     & diagp,diagq,diagr,
     & dimdiagp,dimdiagq,dimdiagr,addp,addq,addr)
c
c     this routine realize following procedure
c     for symp.ne.symq.ne.symr
c
c     ec = sum(p,q,r) [ W(p,q,r) . V(p,q,r) / Dijkpqr ]
c
c     w      - W  matrix (I)
c     v      - V  matrix (I)
c     dimp   - dimension of p index (I)
c     dimq   - dimension of q index (I)
c     dimr   - dimension of r index (I)
c     denijk - sum of i,j,k diagonal (other) F parts (I)
c     ec     - energy contribution
c     diagp  - vector of diagonal parts of symmetry p (I)
c     diagq  - vector of diagonal parts of symmetry q (I)
c     diagr  - vector of diagonal parts of symmetry r (I)
c     dimdiagp - dimension of diagonal parts in symetry p (I)
c     dimdiagq - dimension of diagonal parts in symetry q (I)
c     dimdiagr - dimension of diagonal parts in symetry r (I)
c     addp   - additional constant to p (Nocc) (I)
c     addq   - additional constant to q (Nocc) (I)
c     addr   - additional constant to r (Nocc) (I)
c
       integer dimp,dimq,dimr
       integer dimdiagp,dimdiagq,dimdiagr,addp,addq,addr
       real*8 ec,denijk
       real*8 w(1:dimp,1:dimq,1:dimr)
       real*8 v(1:dimp,1:dimq,1:dimr)
       real*8 diagp(1:dimdiagp)
       real*8 diagq(1:dimdiagq)
       real*8 diagr(1:dimdiagr+dimr)
c
c     help variables
c
       integer p,q,r
       real*8 denijkr,denijkqr,denijkpqr
c
       ec=0.0d0
c
       do 100 r=1,dimr
       denijkr=denijk-diagr(addr+r)
       do 101 q=1,dimq
       denijkqr=denijkr-diagq(q+addq)
       do 102 p=1,dimp
       denijkpqr=denijkqr-diagp(p+addp)
       ec=ec+w(p,q,r)*v(p,q,r)/denijkpqr
 102    continue
 101    continue
 100    continue
c
       return
       end
