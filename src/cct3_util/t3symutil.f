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
c
c     this file contains following routines:
c
c     t3grc0
c     minusa
c     setb
c     stz
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
       subroutine cct3_t3grc0 (nind,typ,typp,typq,typr,typs,stot,
     & poss0,posst,mapd,mapi)
c
c     nind   - number of indexes (I)
c     typ    - typ of mediate (I)
c     typp   - typ of index p (I)
c     typq   - typ of index q (I)
c     typr   - typ of index r (I)
c     typs   - typ of index s (I)
c     stot   - overall symetry of the mediate (I)
c     poss0  - initil possition of mediate (I)
c     posst  - final possition of the mediate (O)
c     mapd   - direct map of the mediate (O)
c     mapi   - inverse map of the mediate (O)
c
c     this routine defines mapd and mapi for given intermediat
c     it can done exactly the same maps like grc0 in CCSD
c     plus additional types of mediates are introduced:
c     type    meaning
c     5     p>q>r,s ; also p>q>r
c     6     p,q>r>s
c     7     p>=q,r,s ; also p>=q,r; p>=q
c     8          p,q>=r,s ; also p,q>=s
c     9     p,q,q>=s
c     10     p>=q,r>=s
c     11     p>=q>=r,s ; also p>=q>=r
c     12     p,q>=r>=s
c
c     currently, these new types are implemented only for nind=3
c
c     !N.B. (this routine cannot run with +OP2)
c     N.B. this routine do not test stupidities
c
c
       integer nind,typ,typp,typq,typr,typs,stot,poss0,posst
c
#include "t31.fh"
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
c     help variables
c
       integer sp,sq,sr,ss,spq,spqr
       integer nsymq,nsymr
       integer poss,i,nhelp1,nhelp2,nhelp3,nhelp4
       integer rsk1,rsk2

c     To fix some compiler warnings

      ss=0
      poss=0
      rsk1=0
      rsk2=0
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
c     def reucion sumations keys : rsk1 for pq, rsk2 for qr
c
       if (typ.eq.0) then
       rsk1=0
       rsk2=0
       else if (typ.eq.1) then
       rsk1=1
       rsk2=0
       else if (typ.eq.2) then
       rsk1=0
       rsk2=1
       else if (typ.eq.5) then
       rsk1=1
       rsk2=1
       else if (typ.eq.7) then
       rsk1=1
       rsk2=0
       else if (typ.eq.8) then
       rsk1=0
       rsk2=1
       else if (typ.eq.11) then
       rsk1=1
       rsk2=1
       end if
c
       i=1
       poss=poss0
c
       do 200 sp=1,nsym
       if (rsk1.eq.1) then
       nsymq=sp
       else
       nsymq=nsym
       end if
c
       do 201 sq=1,nsymq
       spq=mmul(sp,sq)
c
       sr=mmul(stot,spq)
       if ((rsk2.eq.1).and.(sq.lt.sr)) then
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
       else if (typ.eq.5) then
       if (sp.eq.sr) then
       mapd(i,2)=nhelp1*(nhelp1-1)*(nhelp1-2)/6
       else if (sp.eq.sq) then
       mapd(i,2)=nhelp1*(nhelp1-1)*nhelp3/2
       else if (sq.eq.sr) then
       mapd(i,2)=nhelp1*nhelp2*(nhelp2-1)/2
       else
       mapd(i,2)=nhelp1*nhelp2*nhelp3
       end if
       else if ((typ.eq.7).and.(sp.eq.sq)) then
       mapd(i,2)=nhelp1*(nhelp1+1)*nhelp3/2
       else if ((typ.eq.8).and.(sq.eq.sr)) then
       mapd(i,2)=nhelp1*nhelp2*(nhelp2+1)/2
       else if (typ.eq.11) then
       if (sp.eq.ss) then
       mapd(i,2)=nhelp1*(nhelp1+1)*(nhelp1+2)/6
       else if (sp.eq.sq) then
       mapd(i,2)=nhelp1*(nhelp1+1)*nhelp3/2
       else if (sq.eq.sr) then
       mapd(i,2)=nhelp1*nhelp2*(nhelp2+1)/2
       else
       mapd(i,2)=nhelp1*nhelp2*nhelp3
       end if
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
c
c     -----------------------------
c
       subroutine minusa (wrk,wrksize,
     & mapda,factor)
c
c     this routine do
c     A = factor . A
c
c     mapda  - direct map of A m(I/O)
c     factor - numerical factor (I)
c
c     N.B. this routine should be done using matrix operations

#include "wrk.fh"
       integer mapda(0:512,1:6)
       real*8 factor
c
c     help variables
c
       integer nhelp1,nhelp2,nhelp3
c
c
c1    def the length of the mediate
       nhelp1=mapda(0,5)
       nhelp3=mapda(nhelp1,1)+mapda(nhelp1,2)-mapda(1,1)
c
c2    def initial possition
       nhelp2=mapda(1,1)
c
c3    refactoring
       do 100 nhelp1=nhelp2,nhelp2+nhelp3-1
       wrk(nhelp1)=factor*wrk(nhelp1)
 100    continue
c
       return
       end
c
c     -----------------------------
c
       subroutine setb (wrk,wrksize,
     & mapda,mapdb,factor)
c
c     this routine do
c     B = factor . A
c
c     mapda  - direct map of A (I)
c     mapdb  - direct map of B (I)
c     factor - numerical factor (I)
c
c     mediate B must have defined maps, and their must be
c     of identicat type as those for A. If they are not defined,
c     use grc0 before setb
c
c     N.B. this routine should be done using matrix operations

#include "wrk.fh"
       integer mapda(0:512,1:6)
       integer mapdb(0:512,1:6)
       real*8 factor
c
c     help variables
c
       integer possa0,possb0,length,nhelp
c
c
c1    def the length of the mediate
       nhelp=mapda(0,5)
       length=mapda(nhelp,1)+mapda(nhelp,2)-mapda(1,1)
       if (length.eq.0) return
c
c2    def initial possitions
       possa0=mapda(1,1)
       possb0=mapdb(1,1)
c
c3    set B=f.A
       do 100 nhelp=0,length-1
       wrk(possb0+nhelp)=factor*wrk(possa0+nhelp)
 100    continue
c
       return
       end
c
c     -------------------------------------
c
       subroutine stz (wrk,wrksize,
     & mapda)
c
c     this routine vanish A
c     A = 0
c
c     mapda  - direct map of A m(I/O)
c
c     N.B. this routine should be done using matrix operations

#include "wrk.fh"
       integer mapda(0:512,1:6)
c
c     help variables
c
       integer nhelp1,nhelp2,nhelp3
c
c
c1    def the length of the mediate
       nhelp1=mapda(0,5)
       nhelp3=mapda(nhelp1,1)+mapda(nhelp1,2)-mapda(1,1)
c
c2    def initial possition
       nhelp2=mapda(1,1)
c
c3    refactoring
       do 100 nhelp1=nhelp2,nhelp2+nhelp3-1
       wrk(nhelp1)=0.0d0
 100    continue
c
       return
       end
c
c     -------------------------------------
c
