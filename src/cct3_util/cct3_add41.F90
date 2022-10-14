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
       subroutine cct3_add41 (a,b,p,dimp,dimpq,dimr,fact)

!     this routine do:
!     B(pq,r) <-- fact * A(q,r) for given p
!
#include "t31.fh"
       integer dimp,dimpq,dimr,p
       real*8 fact
       real*8 b(1:dimpq,1:dimr)
       real*8 a(1:dimp,1:dimr)
!
!     help variable
!
       integer q,r,pq,qp
!
       if (p.eq.1) goto 101
!
       do 100 r=1,dimr
       pq=nshf(p)
!
       do 50 q=1,p-1
       pq=pq+1
       b(pq,r)=b(pq,r)+fact*a(q,r)
 50     continue
!
 100    continue
!
 101    if (p.eq.dimp) then
       return
       end if
!
       do 200 r=1,dimr
!
       do 150 q=p+1,dimp
       qp=nshf(q)+p
       b(qp,r)=b(qp,r)-fact*a(q,r)
 150    continue
!
 200    continue
!
       return
       end
