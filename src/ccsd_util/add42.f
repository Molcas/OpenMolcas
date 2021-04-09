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
       subroutine add42 (a,b,q,dimq,dimpq,dimr,fact)

c     this routine do:
c     B(pq,r) <-- fact * A(p,r) for given q
c
#include "ccsd1.fh"
       integer dimq,dimpq,dimr,q
       real*8 fact
       real*8 b(1:dimpq,1:dimr)
       real*8 a(1:dimq,1:dimr)
c
c     help variable
c
       integer pq,qp,r,p
c
       if (q.eq.1) goto 101
c
       do 100 r=1,dimr
       qp=nshf(q)
c
       do 50 p=1,q-1
       qp=qp+1
       b(qp,r)=b(qp,r)-fact*a(p,r)
 50     continue
c
 100    continue
c
 101    if (q.eq.dimq) then
       return
       end if
c
       do 200 r=1,dimr
c
       do 150 p=q+1,dimq
       pq=nshf(p)+q
       b(pq,r)=b(pq,r)+fact*a(p,r)
 150    continue
c
 200    continue
c
       return
       end
