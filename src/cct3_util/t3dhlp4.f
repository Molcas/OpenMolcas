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
       subroutine t3dhlp4 (w,v,dimp,dimpqr,denijk,ec,                   &
     & diagp,                                                           &
     & dimdiagp,addp)
!
!     this routine realize following procedure
!     for symp=symq=symr
!
!     ec = sum(pqr) [ W(pqr) . V(pqr) / Dijkpqr ]
!
!     w      - W  matrix (I)
!     v      - V  matrix (I)
!     dimp   - dimension of p index (I)
!     dimpqr - dimension of q index (I)
!     denijk - sum of i,j,k diagonal (other) F parts (I)
!     ec     - energy contribution
!     diagp  - vector of diagonal parts of symmetry p (I)
!     dimdiagp - dimension of diagonal parts in symetry p (I)
!     addp   - additional constant to p (Nocc) (I)
!
       integer dimp,dimpqr
       integer dimdiagp,addp
       real*8 ec,denijk
       real*8 w(1:dimpqr)
       real*8 v(1:dimpqr)
       real*8 diagp(1:dimdiagp)
!
!     help variables
!
       integer p,q,r,pqr
       real*8 denijkp,denijkpq,denijkpqr
!
       ec=0.0d0
!
       pqr=0
       do 100 p=3,dimp
       denijkp=denijk-diagp(p+addp)
       do 101 q=2,p-1
       denijkpq=denijkp-diagp(q+addp)
       do 102 r=1,q-1
       denijkpqr=denijkpq-diagp(r+addp)
       pqr=pqr+1
       ec=ec+w(pqr)*v(pqr)/denijkpqr
 102    continue
 101    continue
 100    continue
!
       return
       end
