!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
       subroutine cct3_map22 (a,b,dimp,dimq,dim1,dim2,p,q,nfact)
!
       integer dimp,dimq,dim1,dim2,p,q,nfact
       real*8 a(1:dimp,1:dimq)
       real*8 b(1:dim1,1:dim2)
!       integer index(1:2)
!
       integer pp,qq
!
       if (nfact.eq.1) then
!
!     nfact = + 1
!
!     do 100 qq=1,dimq
!     index(q)=qq
!     do 100 pp=1,dimp
!     index(p)=pp
!     b(index(1),index(2))=a(pp,qq)
!100  continue
!
       if (p.eq.1) then
!     1* (2)
       do 110 qq=1,dimq
       do 111 pp=1,dimp
       b(pp,qq)=a(pp,qq)
 111    continue
 110    continue
       else
!     2* (1)
       do 120 qq=1,dimq
       do 121 pp=1,dimp
       b(qq,pp)=a(pp,qq)
 121    continue
 120    continue
       end if
!
!
       else
!
!     nfact = -1
!
!     do 200 qq=1,dimq
!     index(q)=qq
!     do 200 pp=1,dimp
!     index(p)=pp
!     b(index(1),index(2))=-a(pp,qq)
!200  continue
!
       if (p.eq.1) then
!     1* (2)
       do 210 qq=1,dimq
       do 211 pp=1,dimp
       b(pp,qq)=-a(pp,qq)
 211    continue
 210    continue
       else
!     2* (1)
       do 220 qq=1,dimq
       do 221 pp=1,dimp
       b(qq,pp)=-a(pp,qq)
 221    continue
 220    continue
       end if
!
       end if
!
       return
! Avoid unused argument warnings
      if (.false.) call Unused_integer(q)
       end
