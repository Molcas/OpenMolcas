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
       subroutine pack211 (ap,am,b,dimp,dimq,rc)
c
c     this routine do: B(p,q) = A+(p,q) - A-(q,p) for symp>symq
c
       integer dimp,dimq,rc
       real*8 ap(1:dimp,1:dimq)
       real*8 am(1:dimq,1:dimp)
       real*8 b(1:dimp,1:dimq)
c
c     help variables
c
       integer p,q
c
       rc=0
       do 100 q=1,dimq
       do 101 p=1,dimp
       b(p,q)=ap(p,q)-am(q,p)
 101    continue
 100    continue
c
       return
       end
