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
       subroutine add44 (a,b,r,dimp,dimqr,dimq,fact)

c     this routine do:
c     B(p,qr) <-- fact * A(p,q) for given r
c
#include "ccsd1.fh"
       integer dimp,dimqr,dimq,r
       real*8 fact
       real*8 b(1:dimp,1:dimqr)
       real*8 a(1:dimp,1:dimq)
c
c     help variable
c
       integer p,qr,rq,q
c
       if (r.eq.1) goto 101
c
       rq=nshf(r)
       do 100 q=1,r-1
       rq=rq+1
c
       do 50 p=1,dimp
       b(p,rq)=b(p,rq)-fact*a(p,q)
 50     continue
c
 100    continue
c
 101    if (r.eq.dimq) then
       return
       end if
c
c
       do 200 q=r+1,dimq
       qr=nshf(q)+r
       do 150 p=1,dimp
       b(p,qr)=b(p,qr)+fact*a(p,q)
 150    continue
c
 200    continue
c
       return
       end
