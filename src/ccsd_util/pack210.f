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
       subroutine pack210 (a,b,dimpq,dimp,rc)
c
c     this routine do: B(pq) = A(p,q) - A(q,p) for symp=symq
c
       integer dimp,dimpq,rc
       real*8 a(1:dimp,1:dimp)
       real*8 b(1:dimpq)
c
c     help variables
c
       integer p,q,pq
c
       rc=0
       if (dimp.gt.1) then
c
       pq=0
       do 100 p=2,dimp
       do 101 q=1,p-1
       pq=pq+1
       b(pq)=a(p,q)-a(q,p)
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
