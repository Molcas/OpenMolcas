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
       subroutine expand0 (a,b,dimpq,dimp)
c
c     expand a(pq) -> b(p,q)
c     assumption: p>q, a(p,q)=-a(q,p)
c     RISC version
c
       integer dimpq,dimp
       real*8 a(1:dimpq)
       real*8 b(1:dimp,1:dimp)
c
c     help variables
c
       integer p,q,pq
       real*8 scalar
c
       if (dimp.gt.1) then
c
       pq=0
       do 100 p=2,dimp
       do 101 q=1,p-1
       pq=pq+1
       scalar=a(pq)
       b(p,q)=scalar
       b(q,p)=-scalar
 101    continue
 100    continue
c
       end if
c
       do p=1,dimp
       b(p,p)=0.0d0
       end do
c
       return
       end
