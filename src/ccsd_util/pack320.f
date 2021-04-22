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
       subroutine pack320 (a,b,dimp,dimqr,dimq,rc)
c
c     this routine do: B(p,qr) = A(p,q,r) - A(p,r,q) for symq=symr
c
       integer dimp,dimqr,dimq,rc
       real*8 a(1:dimp,1:dimq,1:dimq)
       real*8 b(1:dimp,1:dimqr)
c
c     help variables
c
       integer p,q,r,qr
c
       rc=0
       if (dimq.gt.1) then
c
       qr=0
       do 100 q=2,dimq
       do 101 r=1,q-1
       qr=qr+1
       do 102 p=1,dimp
       b(p,qr)=a(p,q,r)-a(p,r,q)
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
