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
       subroutine map41 (a,b,dimp,dimq,dimr,dims,p,q,r,s,nfact)
c
c     maping A(p1,q1,r1,s1) -> nfact* B(p2,q2,r2,s2)
c
       real*8 a(*)
       real*8 b(*)
       integer dim(4)
       integer dimp,dimq,dimr,dims,p,q,r,s,nfact
       dim(p)=dimp
       dim(q)=dimq
       dim(r)=dimr
       dim(s)=dims
       call map42 (a,b,dimp,dimq,dimr,dims,dim(1),dim(2),dim(3),dim(4),
     &             p,q,r,s,nfact)
       return
       end
