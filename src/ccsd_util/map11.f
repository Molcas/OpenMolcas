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
       subroutine map11 (a,b,dimp,nfact)
c
c     maping A(p) -> nfact*B(p)
c
       integer dimp,nfact
       real*8 a(1:dimp)
       real*8 b(1:dimp)
c
       integer pp
c
       if (nfact.eq.1) then
c
       do 100 pp=1,dimp
       b(pp)=a(pp)
 100    continue
c
       else
c
       do 200 pp=1,dimp
       b(pp)=-a(pp)
 200    continue
c
       end if
c
       return
       end
