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
       subroutine pack311 (ap,am,b,dimp,dimq,dimr,rc)
c
c     this routine do: B(p,q,r) = A+(p,q,r) - A-(q,p,r) for symp>symq
c
       integer dimp,dimq,dimr,rc
       real*8 ap(1:dimp,1:dimq,1:dimr)
       real*8 am(1:dimq,1:dimp,1:dimr)
       real*8 b(1:dimp,1:dimq,1:dimr)
c
c     help variables
c
       integer p,q,r
c
       rc=0
       do 100 r=1,dimr
       do 101 q=1,dimq
       do 102 p=1,dimp
       b(p,q,r)=ap(p,q,r)-am(q,p,r)
 102    continue
 101    continue
 100    continue
c
       return
       end
