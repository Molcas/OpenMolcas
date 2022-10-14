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
       subroutine cct3_add42 (a,b,q,dimq,dimpq,dimr,fact)

!     this routine do:
!     B(pq,r) <-- fact * A(p,r) for given q
!
#include "t31.fh"
       integer dimq,dimpq,dimr,q
       real*8 fact
       real*8 b(1:dimpq,1:dimr)
       real*8 a(1:dimq,1:dimr)
!
!     help variable
!
       integer pq,qp,r,p
!
       if (q.eq.1) goto 101
!
       do 100 r=1,dimr
       qp=nshf(q)
!
       do 50 p=1,q-1
       qp=qp+1
       b(qp,r)=b(qp,r)-fact*a(p,r)
 50     continue
!
 100    continue
!
 101    if (q.eq.dimq) then
       return
       end if
!
       do 200 r=1,dimr
!
       do 150 p=q+1,dimq
       pq=nshf(p)+q
       b(pq,r)=b(pq,r)+fact*a(p,r)
 150    continue
!
 200    continue
!
       return
       end
