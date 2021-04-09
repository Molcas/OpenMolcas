************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2006, Pavel Neogrady                                   *
************************************************************************
       subroutine grc0stack (bsize,typ,typp,typq,typr,typs,stot,
     & poss0,posst,mapd,mapi)
c
c             This routine defines mapd and mapi for specific
c        3 index intermediat A(pq,Bp), needed when stacking
c        (About Bp, see notes in multstack)
c        This routine is a modification of grc0 routine
c
c        P.N. 17.02.06
c     !N.B. (this routine cannot run with +OP2)
c
       integer bsize,typ,typp,typq,typr,typs,stot,poss0,posst
c
#include "ccsd1.fh"
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
c     help variables
c
       integer sp,sq
CLD    integer sp,sq,sr,ss,spq,spqr
CLD    integer nsymq,nsymr
       integer poss,i,nhelp1,nhelp2,nhelp3
CLD    integer poss,i,nhelp1,nhelp2,nhelp3,nhelp4

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
c     matrix A(p,q) or specifilally A(i,j,Bp)
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
       mapd(i,2)=bsize*nhelp1*(nhelp1-1)/2
       else
       mapd(i,2)=bsize*nhelp1*nhelp2
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
