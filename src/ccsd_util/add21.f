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
       subroutine add21 (a,b,p,dimp,dimq,fact)
c
c     this routine do:
c     B(p,q) <-- fact * A(q) for given p
c
       integer dimp,dimq,p
       real*8 fact
       real*8 b(1:dimp,1:dimq)
       real*8 a(1:dimq)
c
c     help variable
c
       integer q
c
       do 100 q=1,dimq
       b(p,q)=b(p,q)+fact*a(q)
 100    continue
c
       return
       end
