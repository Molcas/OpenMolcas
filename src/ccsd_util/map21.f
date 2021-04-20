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
       subroutine map21 (a,b,dimp,dimq,p,q,nfact)
c
c     maping A(p1,q1) -> nfact*B(p2,q2)
c
       real*8 a(*)
       real*8 b(*)
       integer dim(2)
       integer dimp,dimq,p,q,nfact
       dim(p)=dimp
       dim(q)=dimq
       call map22 (a,b,dimp,dimq,dim(1),dim(2),p,q,nfact)
       return
       end
