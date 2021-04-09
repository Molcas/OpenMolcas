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
       subroutine t3dhlp3 (w,v,dimp,dimq,dimqr,denijk,ec,
     & diagp,diagq,
     & dimdiagp,dimdiagq,addp,addq)
c
c     this routine realize following procedure
c     for symp.ne.symq=symr
c
c     ec = sum(p,qr) [ W(p,qr) . V(p,qr) / Dijkpqr ]
c
c     w      - W  matrix (I)
c     v      - V  matrix (I)
c     dimp   - dimension of p index (I)
c     dimq   - dimension of q (r) index (I)
c     dimqr  - dimension of r index (I)
c     denijk - sum of i,j,k diagonal (other) F parts (I)
c     ec     - energy contribution
c     diagp  - vector of diagonal parts of symmetry p (I)
c     diagq  - vector of diagonal parts of symmetry q (I)
c     dimdiagp - dimension of diagonal parts in symetry p (I)
c     dimdiagq - dimension of diagonal parts in symetry q (I)
c     addp   - additional constant to p (Nocc) (I)
c     addq   - additional constant to q (Nocc) (I)
c
       integer dimp,dimq,dimqr
       integer dimdiagp,dimdiagq,addp,addq
       real*8 ec,denijk
       real*8 w(1:dimp,1:dimqr)
       real*8 v(1:dimp,1:dimqr)
       real*8 diagp(1:dimdiagp)
       real*8 diagq(1:dimdiagq)
c
c     help variables
c
       integer p,q,r,qr
       real*8 denijkq,denijkqr,denijkpqr
c
       ec=0.0d0
c
       qr=0
       do 100 q=2,dimq
       denijkq=denijk-diagq(addq+q)
       do 101 r=1,q-1
       denijkqr=denijkq-diagq(r+addq)
       qr=qr+1
       do 102 p=1,dimp
       denijkpqr=denijkqr-diagp(p+addp)
       ec=ec+w(p,qr)*v(p,qr)/denijkpqr
 102    continue
 101    continue
 100    continue
c
       return
       end
