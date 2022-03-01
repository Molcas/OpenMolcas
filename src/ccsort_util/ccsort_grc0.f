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
       subroutine ccsort_grc0 (nind,typ,typp,typq,typr,typs,stot,       &
     & poss0,posst,mapd,mapi)
!
!     this routine defines mapd and mapi for given intermediat
!
#include "ccsort.fh"
       integer nind,typ,typp,typq,typr,typs,stot,poss0,posst
!
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
!     integer mul(1:8,1:8)
       integer dimm(1:4,1:8)
!
!     help variables
!
       integer sp,sq,sr,ss,spq,spqr
       integer nsymq,nsymr
       integer poss,i,nhelp1,nhelp2,nhelp3,nhelp4
!
!*    !!!!!!!! def dimm to je tu len terazky !!!!
!
       do i=1,nsym
       dimm(1,i)=noa(i)
       dimm(2,i)=nob(i)
       dimm(3,i)=nva(i)
       dimm(4,i)=nvb(i)
       end do
!
!     vanishing mapi files
!
       do nhelp1=1,nsym
       do nhelp2=1,nsym
       do nhelp3=1,nsym
       mapi(nhelp3,nhelp2,nhelp1)=0
       end do
       end do
       end do
!
!      To get rid of tiring compiler warning
       poss=0
       if (nind.eq.1) then
!
!     matrix A(p)
!
       i=1
       poss=poss0
       sp=mul(stot,1)
!
       nhelp1=dimm(typp,sp)
!
!     def mapi
       mapi(1,1,1)=i
!
!     def possition
       mapd(i,1)=poss
!
!     def length
       mapd(i,2)=nhelp1
!
!     def sym p,q
       mapd(i,3)=sp
       mapd(i,4)=0
       mapd(i,5)=0
       mapd(i,6)=0
!
       poss=poss+mapd(i,2)
       i=i+1
!
       else if (nind.eq.2) then
!
!     matrix A(p,q)
!
       i=1
       poss=poss0
!
       do 100 sp=1,nsym
!
       sq=mul(stot,sp)
       if ((typ.eq.1).and.(sp.lt.sq)) then
!     Meggie out
       goto 100
       end if
!
       nhelp1=dimm(typp,sp)
       nhelp2=dimm(typq,sq)
!
!     def mapi
       mapi(sp,1,1)=i
!
!     def possition
       mapd(i,1)=poss
!
!     def length
       if ((typ.eq.1).and.(sp.eq.sq)) then
       mapd(i,2)=nhelp1*(nhelp1-1)/2
       else
       mapd(i,2)=nhelp1*nhelp2
       end if
!
!     def sym p,q
       mapd(i,3)=sp
       mapd(i,4)=sq
       mapd(i,5)=0
       mapd(i,6)=0
!
       poss=poss+mapd(i,2)
       i=i+1
!
 100    continue
!
       else if (nind.eq.3) then
!
!     matrix A(p,q,r)
!
       i=1
       poss=poss0
!
       do 200 sp=1,nsym
       if (typ.eq.1) then
       nsymq=sp
       else
       nsymq=nsym
       end if
!
       do 201 sq=1,nsymq
       spq=mul(sp,sq)
!
       sr=mul(stot,spq)
       if ((typ.eq.2).and.(sq.lt.sr)) then
!     Meggie out
       goto 201
       end if
!
       nhelp1=dimm(typp,sp)
       nhelp2=dimm(typq,sq)
       nhelp3=dimm(typr,sr)
!
!     def mapi
       mapi(sp,sq,1)=i
!
!     def possition
       mapd(i,1)=poss
!
!     def length
       if ((typ.eq.1).and.(sp.eq.sq)) then
       mapd(i,2)=nhelp1*(nhelp1-1)*nhelp3/2
       else if ((typ.eq.2).and.(sq.eq.sr)) then
       mapd(i,2)=nhelp1*nhelp2*(nhelp2-1)/2
       else
       mapd(i,2)=nhelp1*nhelp2*nhelp3
       end if
!
!     def sym p,q,r
       mapd(i,3)=sp
       mapd(i,4)=sq
       mapd(i,5)=sr
       mapd(i,6)=0
!
       poss=poss+mapd(i,2)
       i=i+1
!
 201    continue
 200    continue
!
       else if (nind.eq.4) then
!
!     matrix A(p,q,r,s)
!
       i=1
       poss=poss0
!
       do 300 sp=1,nsym
       if ((typ.eq.1).or.(typ.eq.4)) then
       nsymq=sp
       else
       nsymq=nsym
       end if
!
       do 301 sq=1,nsymq
       spq=mul(sp,sq)
       if (typ.eq.2) then
       nsymr=sq
       else
       nsymr=nsym
       end if
!
       do 302 sr=1,nsymr
       spqr=mul(spq,sr)
!
       ss=mul(stot,spqr)
       if (((typ.eq.3).or.(typ.eq.4)).and.(sr.lt.ss)) then
!     Meggie out
       goto 302
       end if
!
       nhelp1=dimm(typp,sp)
       nhelp2=dimm(typq,sq)
       nhelp3=dimm(typr,sr)
       nhelp4=dimm(typs,ss)
!
!     def mapi
       mapi(sp,sq,sr)=i
!
!     def possition
       mapd(i,1)=poss
!
!     def length
       if ((typ.eq.1).and.(sp.eq.sq)) then
       mapd(i,2)=nhelp1*(nhelp2-1)*nhelp3*nhelp4/2
       else if ((typ.eq.2).and.(sq.eq.sr)) then
       mapd(i,2)=nhelp1*nhelp2*(nhelp3-1)*nhelp4/2
       else if ((typ.eq.3).and.(sr.eq.ss)) then
       mapd(i,2)=nhelp1*nhelp2*nhelp3*(nhelp4-1)/2
       else if (typ.eq.4) then
       if ((sp.eq.sq).and.(sr.eq.ss)) then
       mapd(i,2)=nhelp1*(nhelp2-1)*nhelp3*(nhelp4-1)/4
       else if (sp.eq.sq) then
       mapd(i,2)=nhelp1*(nhelp2-1)*nhelp3*nhelp4/2
       else if (sr.eq.ss) then
       mapd(i,2)=nhelp1*nhelp2*nhelp3*(nhelp4-1)/2
       else
       mapd(i,2)=nhelp1*nhelp2*nhelp3*nhelp4
       end if
       else
       mapd(i,2)=nhelp1*nhelp2*nhelp3*nhelp4
       end if
!
!     def sym p,q,r,s
       mapd(i,3)=sp
       mapd(i,4)=sq
       mapd(i,5)=sr
       mapd(i,6)=ss
!
       poss=poss+mapd(i,2)
       i=i+1
!
 302    continue
 301    continue
 300    continue
!
       end if

!
       posst=poss
!
!     definition of other coll
!
       mapd(0,1)=typp
       mapd(0,2)=typq
       mapd(0,3)=typr
       mapd(0,4)=typs
       mapd(0,5)=i-1
       mapd(0,6)=typ
!
       return
       end
