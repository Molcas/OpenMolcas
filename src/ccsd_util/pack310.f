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
       subroutine pack310 (a,b,dimpq,dimr,dimp,rc)
c
c     this routine do: B(pq,r) = A(p,q,r) - A(q,p,r) for symp=symq
c
       integer dimp,dimpq,dimr,rc
       real*8 a(1:dimp,1:dimp,1:dimr)
       real*8 b(1:dimpq,1:dimr)
c
c     help variables
c
       integer p,q,r,pq
c
       rc=0
       if (dimp.gt.1) then
c
       do 100 r=1,dimr
       pq=0
       do 101 p=2,dimp
       do 102 q=1,p-1
       pq=pq+1
       b(pq,r)=a(p,q,r)-a(q,p,r)
 102    continue
 101    continue
 100    continue
c
       else
c     RC=1 : dimp is less than 2
       rc=1
       end if
c
       return
       end
