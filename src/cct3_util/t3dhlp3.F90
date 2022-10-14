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
       subroutine t3dhlp3 (w,v,dimp,dimq,dimqr,denijk,ec,               &
     & diagp,diagq,                                                     &
     & dimdiagp,dimdiagq,addp,addq)
!
!     this routine realize following procedure
!     for symp.ne.symq=symr
!
!     ec = sum(p,qr) [ W(p,qr) . V(p,qr) / Dijkpqr ]
!
!     w      - W  matrix (I)
!     v      - V  matrix (I)
!     dimp   - dimension of p index (I)
!     dimq   - dimension of q (r) index (I)
!     dimqr  - dimension of r index (I)
!     denijk - sum of i,j,k diagonal (other) F parts (I)
!     ec     - energy contribution
!     diagp  - vector of diagonal parts of symmetry p (I)
!     diagq  - vector of diagonal parts of symmetry q (I)
!     dimdiagp - dimension of diagonal parts in symetry p (I)
!     dimdiagq - dimension of diagonal parts in symetry q (I)
!     addp   - additional constant to p (Nocc) (I)
!     addq   - additional constant to q (Nocc) (I)
!
       integer dimp,dimq,dimqr
       integer dimdiagp,dimdiagq,addp,addq
       real*8 ec,denijk
       real*8 w(1:dimp,1:dimqr)
       real*8 v(1:dimp,1:dimqr)
       real*8 diagp(1:dimdiagp)
       real*8 diagq(1:dimdiagq)
!
!     help variables
!
       integer p,q,r,qr
       real*8 denijkq,denijkqr,denijkpqr
!
       ec=0.0d0
!
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
!
       return
       end
