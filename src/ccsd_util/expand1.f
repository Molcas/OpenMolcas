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
       subroutine expand1 (a,b,dimpq,dimr,dimp)
c
c     expand a(pq,r) -> b(p,q,r)
c     assumption: p>q, a(p,q,r)=-a(q,p,r)
c     RISC version
c
       integer dimpq,dimr,dimp
       real*8 a(1:dimpq,1:dimr)
       real*8 b(1:dimp,1:dimp,1:dimr)
c
c     help variables
c
       integer p,q,r,pq
       real*8 scalar
c
       if (dimp.gt.1) then
c
       do r=1,dimr
       pq=0
       do 100 p=2,dimp
       do 101 q=1,p-1
       pq=pq+1
       scalar=a(pq,r)
       b(p,q,r)=scalar
       b(q,p,r)=-scalar
 101    continue
 100    continue
       end do
c
       end if
c
       do r=1,dimr
       do p=1,dimp
       b(p,p,r)=0.0d0
       end do
       end do
c
       return
       end
