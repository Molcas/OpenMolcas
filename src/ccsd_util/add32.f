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
       subroutine add32 (a,b,q,dimp,dimq,dimr,fact)
c
c     this routine do:
c     B(p,q,r) <-- fact * A(p,r) for given q
c
       integer dimp,dimq,dimr,q
       real*8 fact
       real*8 b(1:dimp,1:dimq,1:dimr)
       real*8 a(1:dimp,1:dimr)
c
c     help variable
c
       integer p,r
c
       do 100 r=1,dimr
       do 101 p=1,dimp
       b(p,q,r)=b(p,q,r)+fact*a(p,r)
 101    continue
 100    continue
c
       return
       end
