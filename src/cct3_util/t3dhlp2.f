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
       subroutine t3dhlp2 (w,v,dimp,dimpq,dimr,denijk,ec,
     & diagp,diagr,
     & dimdiagp,dimdiagr,addp,addr)
c
c     this routine realize following procedure
c     for symp=symq.ne.symr
c
c     ec = sum(pq,r) [ W(pq,r) . V(pq,r) / Dijkpqr ]
c
c     w      - W  matrix (I)
c     v      - V  matrix (I)
c     dimp   - dimension of p (q) index (I)
c     dimpq  - dimension of pq index (I)
c     dimr   - dimension of r index (I)
c     denijk - sum of i,j,k diagonal (other) F parts (I)
c     ec     - energy contribution
c     diagp  - vector of diagonal parts of symmetry p (q) (I)
c     diagr  - vector of diagonal parts of symmetry r (I)
c     dimdiagp - dimension of diagonal parts in symetry p (I)
c     dimdiagr - dimension of diagonal parts in symetry r (I)
c     addp   - additional constant to p (Nocc) (I)
c     addr   - additional constant to r (Nocc) (I)
c
       integer dimp,dimpq,dimr
       integer dimdiagp,dimdiagr,addp,addr
       real*8 ec,denijk
       real*8 w(1:dimpq,1:dimr)
       real*8 v(1:dimpq,1:dimr)
       real*8 diagp(1:dimdiagp)
       real*8 diagr(1:dimdiagr)
c
c     help variables
c
       integer p,q,r,pq
       real*8 denijkr,denijkpr,denijkpqr
c
       ec=0.0d0
c
       do 100 r=1,dimr
       denijkr=denijk-diagr(addr+r)
       pq=0
       do 101 p=2,dimp
       denijkpr=denijkr-diagp(p+addp)
       do 102 q=1,p-1
       pq=pq+1
       denijkpqr=denijkpr-diagp(q+addp)
       ec=ec+w(pq,r)*v(pq,r)/denijkpqr
 102    continue
 101    continue
 100    continue
c
       return
       end
