!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994, Per Ake Malmqvist                                *
!               1995,1996, Pavel Neogrady                              *
!***********************************************************************
       subroutine mod1 (nsym,nfro,nish,nash,nssh,ndel,norb,nfror,ndelr, &
     & firas,fi,epsras,eps)
!
!     this routine do:
!     1) reduce firas, epsras if nfror>nfro, ndelr>ndel
!     2) redefine nfro,nish,nash,nssh,ndel,norb to proper ones
!
       integer nsym
       integer nfro(1:8)
       integer ndel(1:8)
       integer nish(1:8)
       integer nash(1:8)
       integer nssh(1:8)
       integer norb(1:8)
       integer nfror(1:8)
       integer ndelr(1:8)
       real*8 fi(*)
       real*8 firas(*)
       real*8 eps(*)
       real*8 epsras(*)
!
!     help variables
!
       integer p,q,pqras,pqnew,pras,pnew,isym
       integer ndf,ndd,nup,nlow

!
!1    reduce fi
!
       pqras=0
       pqnew=0
       do 100 isym=1,nsym
!
       ndf=nfror(isym)-nfro(isym)
       ndd=ndelr(isym)-ndel(isym)
       nlow=ndf+1
       nup=norb(isym)-ndd
!
       do 50 p=1,norb(isym)
       do 51 q=1,p
       pqras=pqras+1
!
       if ((p.ge.nlow).and.(p.le.nup)) then
       if ((q.ge.nlow).and.(q.le.nup)) then
       pqnew=pqnew+1
       fi(pqnew)=firas(pqras)
       end if
       end if
!
 51     continue
 50     continue
!
 100    continue
!
!2    reduce eps
!
       pras=0
       pnew=0
       do 200 isym=1,nsym
!
       ndf=nfror(isym)-nfro(isym)
       ndd=ndelr(isym)-ndel(isym)
       nlow=ndf+1
       nup=norb(isym)-ndd
!
       do 150 p=1,norb(isym)
       pras=pras+1
!
       if ((p.ge.nlow).and.(p.le.nup)) then
       pnew=pnew+1
       eps(pnew)=epsras(pras)
       end if
!
 150    continue
!
 200    continue
!
!3    define new nfro,nish,nash,nssh,ndel,norb
!
       do 300 isym=1,nsym
       nash(isym)=nash(isym)
       nish(isym)=nish(isym)-nfror(isym)+nfro(isym)
       nssh(isym)=nssh(isym)-ndelr(isym)+ndel(isym)
       norb(isym)=norb(isym)-nfror(isym)+nfro(isym)-ndelr(isym)         &
     &           +ndel(isym)
       nfro(isym)=nfror(isym)
 300    continue
!
!
       return
       end
