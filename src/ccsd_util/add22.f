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
       subroutine add22 (a,b,q,dimp,dimq,fact)
c
c     this routine do:
c     B(p,q) <-- fact * A(p) for given q
c
       integer dimp,dimq,q
       real*8 fact
       real*8 b(1:dimp,1:dimq)
       real*8 a(1:dimp)
c
c     help variable
c
       integer p
c
       do 100 p=1,dimp
       b(p,q)=b(p,q)+fact*a(p)
 100    continue
c
       return
       end
