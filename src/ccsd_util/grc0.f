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
       subroutine grc0 (nind,typ,typp,typq,typr,typs,stot,
     & poss0,posst,mapd,mapi)
c
c     this routine defines mapd and mapi for given intermediat
c
c     $N.B. (this routine cannot run with +OP2)
c
       integer nind,typ,typp,typq,typr,typs,stot,poss0,posst
c
#include "ccsd1.fh"
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
c     help variables
c
       integer sp,sq,sr,ss,spq,spqr
       integer nsymq,nsymr
       integer poss,i,nhelp1,nhelp2,nhelp3,nhelp4

c     To get rid of compiler warning
      poss=0
      i=0
c
c     vanishing mapi files
c
       do nhelp1=1,nsym
       do nhelp2=1,nsym
       do nhelp3=1,nsym
       mapi(nhelp3,nhelp2,nhelp1)=0
       end do
       end do
       end do
c
       if (nind.eq.1) then
c
c     matrix A(p)
c
       i=1
       poss=poss0
       sp=mmul(stot,1)
c
       nhelp1=dimm(typp,sp)
c
c     def mapi
       mapi(1,1,1)=i
c
c     def possition
       mapd(i,1)=poss
c
c     def length
       mapd(i,2)=nhelp1
c
c     def sym p,q
       mapd(i,3)=sp
       mapd(i,4)=0
       mapd(i,5)=0
       mapd(i,6)=0
c
       poss=poss+mapd(i,2)
       i=i+1
c
       else if (nind.eq.2) then
c
c     matrix A(p,q)
c
       i=1
       poss=poss0
c
       do 100 sp=1,nsym
c
       sq=mmul(stot,sp)
       if ((typ.eq.1).and.(sp.lt.sq)) then
c     Meggie out
       goto 100
       end if
c
       nhelp1=dimm(typp,sp)
       nhelp2=dimm(typq,sq)
c
c     def mapi
       mapi(sp,1,1)=i
c
c     def possition
       mapd(i,1)=poss
c
c     def length
       if ((typ.eq.1).and.(sp.eq.sq)) then
       mapd(i,2)=nhelp1*(nhelp1-1)/2
       else
       mapd(i,2)=nhelp1*nhelp2
       end if
c
c     def sym p,q
       mapd(i,3)=sp
       mapd(i,4)=sq
       mapd(i,5)=0
       mapd(i,6)=0
c
       poss=poss+mapd(i,2)
       i=i+1
c
 100    continue
c
       else if (nind.eq.3) then
c
c     matrix A(p,q,r)
c
       i=1
       poss=poss0
c
       do 200 sp=1,nsym
       if (typ.eq.1) then
       nsymq=sp
       else
       nsymq=nsym
       end if
c
       do 201 sq=1,nsymq
       spq=mmul(sp,sq)
c
       sr=mmul(stot,spq)
       if ((typ.eq.2).and.(sq.lt.sr)) then
c     Meggie out
       goto 201
       end if
c
       nhelp1=dimm(typp,sp)
       nhelp2=dimm(typq,sq)
       nhelp3=dimm(typr,sr)
c
c     def mapi
       mapi(sp,sq,1)=i
c
c     def possition
       mapd(i,1)=poss
c
c     def length
       if ((typ.eq.1).and.(sp.eq.sq)) then
       mapd(i,2)=nhelp1*(nhelp1-1)*nhelp3/2
       else if ((typ.eq.2).and.(sq.eq.sr)) then
       mapd(i,2)=nhelp1*nhelp2*(nhelp2-1)/2
       else
       mapd(i,2)=nhelp1*nhelp2*nhelp3
       end if
c
c     def sym p,q,r
       mapd(i,3)=sp
       mapd(i,4)=sq
       mapd(i,5)=sr
       mapd(i,6)=0
c
       poss=poss+mapd(i,2)
       i=i+1
c
 201    continue
 200    continue
c
       else if (nind.eq.4) then
c
c     matrix A(p,q,r,s)
c
       i=1
       poss=poss0
c
       do 300 sp=1,nsym
       if ((typ.eq.1).or.(typ.eq.4)) then
       nsymq=sp
       else
       nsymq=nsym
       end if
c
       do 301 sq=1,nsymq
       spq=mmul(sp,sq)
       if (typ.eq.2) then
       nsymr=sq
       else
       nsymr=nsym
       end if
c
       do 302 sr=1,nsymr
       spqr=mmul(spq,sr)
c
       ss=mmul(stot,spqr)
       if (((typ.eq.3).or.(typ.eq.4)).and.(sr.lt.ss)) then
c     Meggie out
       goto 302
       end if
c
       nhelp1=dimm(typp,sp)
       nhelp2=dimm(typq,sq)
       nhelp3=dimm(typr,sr)
       nhelp4=dimm(typs,ss)
c
c     def mapi
       mapi(sp,sq,sr)=i
c
c     def possition
       mapd(i,1)=poss
c
c     def length
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
c
c     def sym p,q,r,s
       mapd(i,3)=sp
       mapd(i,4)=sq
       mapd(i,5)=sr
       mapd(i,6)=ss
c
       poss=poss+mapd(i,2)
       i=i+1
c
 302    continue
 301    continue
 300    continue
c
       end if

c
       posst=poss
c
c     definition of other coll
c
       mapd(0,1)=typp
       mapd(0,2)=typq
       mapd(0,3)=typr
       mapd(0,4)=typs
       mapd(0,5)=i-1
       mapd(0,6)=typ
c
       return
       end
