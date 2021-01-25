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
c     grc43y
c     grc42y
c     grc32y
c     grc0
c     multdot
c
c     -------------------------------------------------------
c
       subroutine cct3_grc43y (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,possc0,ix)
c
#include "t31.fh"
c
       integer mapda(0:512,1:6)
       integer mapdb(0:512,1:6)
       integer mapdc(0:512,1:6)
c
       integer mapia(1:8,1:8,1:8)
       integer mapib(1:8,1:8,1:8)
       integer mapic(1:8,1:8,1:8)
c
       integer mvec(1:4096,1:7)
       integer possc0
       integer ssa,ssb
c
c     help variables
c
       integer nhelp1,nhelp2,nhelp4
       integer nhelp41,nhelp42,nhelp43
       integer ntest1,ntest2
       integer sa1,sa2,sa3,sa4,sb1,sb2,sb3,sa23,sa234,sb12
       integer nsymb2
       integer ia,ib,iy,ix
       integer possct
c
c     sctructure A(p,qrs)*B(qrs)=YC(p)
c
c1.0  prepare mapdc,mapic
c
       call cct3_grc0 (1,0,mapda(0,1),0,0,0,mmul(ssa,ssb),
     & possc0,possct,mapdc,mapic)
c
c1.1  define limitations - p,q>r,s must be tested - ntest1
c     p,q,r>s must be tested - ntest2
c
       if (mapda(0,6).eq.2) then
       ntest1=1
       else
       ntest1=0
       end if
c
       if (mapda(0,6).eq.3) then
       ntest2=1
       else
       ntest2=0
       end if
c
c1.2  def symm states and test the limitations
c
       ix=1
       do 100 sb1=1,nsym
       sa2=sb1
       if (ntest1.eq.1) then
       nsymb2=sb1
       else
       nsymb2=nsym
       end if
c
       do 101 sb2=1,nsymb2
       sb12=mmul(sb1,sb2)
       sa3=sb2
       sa23=mmul(sa2,sa3)
c
       sb3=mmul(ssb,sb12)
       sa4=sb3
       sa234=mmul(sa23,sa4)
       if ((ntest2.eq.1).and.(sb2.lt.sb3)) then
c     Meggie out
       goto 101
       end if
c
       sa1=mmul(ssa,sa234)
       if ((ntest1.eq.1).and.(sb1.lt.sb2)) then
c     Meggie out
       goto 101
       end if
c
c1.3  def mvec,mapdc and mapdi
c
       ia=mapia(sa1,sa2,sa3)
       ib=mapib(sb1,sb2,1)
       iy=mapic(1,1,1)
c
c     yes/no
       if ((mapda(ia,2).gt.0).and.(mapdb(ib,2).gt.0)) then
       nhelp1=1
       else
       goto 101
       end if
c
c     rowA
       nhelp2=dimm(mapda(0,1),sa1)
c
c     sum
       nhelp41=dimm(mapda(0,2),sa2)
       nhelp42=dimm(mapda(0,3),sa3)
       nhelp43=dimm(mapda(0,4),sa4)
       if ((ntest1.eq.1).and.(sa2.eq.sa3)) then
       nhelp4=nhelp41*(nhelp41-1)*nhelp43/2
       else if ((ntest2.eq.1).and.(sa3.eq.sa4)) then
       nhelp4=nhelp41*nhelp42*(nhelp42-1)/2
       else
       nhelp4=nhelp41*nhelp42*nhelp43
       end if
c
       mvec(ix,1)=nhelp1
       mvec(ix,2)=mapda(ia,1)
       mvec(ix,3)=mapdb(ib,1)
       mvec(ix,4)=mapdc(iy,1)
       mvec(ix,5)=nhelp2
       mvec(ix,6)=nhelp4
       mvec(ix,7)=0
c
       ix=ix+1
c
 101    continue
 100    continue
       ix=ix-1
c
       return
       end
c
c     --------------------------------
c
       subroutine cct3_grc42y (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,possc0,ix)
c
#include "t31.fh"
c
       integer mapda(0:512,1:6)
       integer mapdb(0:512,1:6)
       integer mapdc(0:512,1:6)
c
       integer mapia(1:8,1:8,1:8)
       integer mapib(1:8,1:8,1:8)
       integer mapic(1:8,1:8,1:8)
c
       integer mvec(1:4096,1:7)
       integer possc0
       integer ssa,ssb
c
c     help variables
c
       integer nhelp1,nhelp2,nhelp4
       integer nhelp41,nhelp42,nhelp21,nhelp22
       integer ntest1,ntest2
       integer sa1,sa2,sa3,sa4,sb1,sb2,sa34,sa134
       integer ia,ib,iy,ix
       integer possct
c
c     sctructure A(pq,rs)*B(rs)=YC(pq)
c
c1.0  prepare mapdc,mapic
c
       if ((mapda(0,6).eq.1).or.(mapda(0,6).eq.4)) then
       ntest1=1
       else
       ntest1=0
       end if
c
       call cct3_grc0 (2,ntest1,mapda(0,1),mapda(0,2),0,0,mmul(ssa,ssb),
     & possc0,possct,mapdc,mapic)
c
c1.1  define limitations - p>q,r,s must be tested - ntest1
c     p,q,r>s must be tested - ntest2
c
       if ((mapda(0,6).eq.1).or.(mapda(0,6).eq.4)) then
       ntest1=1
       else
       ntest1=0
       end if
c
       if ((mapda(0,6).eq.3).or.(mapda(0,6).eq.4)) then
       ntest2=1
       else
       ntest2=0
       end if
c
c1.2  def symm states and test the limitations
c
       ix=1
       do 100 sb1=1,nsym
       sa3=sb1
c
       sb2=mmul(ssb,sb1)
       sa4=sb2
       sa34=mmul(sa3,sa4)
       if ((ntest2.eq.1).and.(sb1.lt.sb2)) then
c     Meggie out
       goto 100
       end if
c
       do 50 sa1=1,nsym
       sa134=mmul(sa1,sa34)
c
       sa2=mmul(ssa,sa134)
       if ((ntest1.eq.1).and.(sa1.lt.sa2)) then
c     Meggie out
       goto 50
       end if
c
c1.3  def mvec,mapdc and mapdi
c
       ia=mapia(sa1,sa2,sa3)
       ib=mapib(sb1,1,1)
       iy=mapic(sa1,1,1)
c
c     yes/no
       if ((mapda(ia,2).gt.0).and.(mapdb(ib,2).gt.0)) then
       nhelp1=1
       else
       goto 50
       end if
c
c     rowA
       nhelp21=dimm(mapda(0,1),sa1)
       nhelp22=dimm(mapda(0,2),sa2)
       if ((ntest1.eq.1).and.(sa1.eq.sa2)) then
       nhelp2=nhelp21*(nhelp21-1)/2
       else
       nhelp2=nhelp21*nhelp22
       end if
c
c     sum
       nhelp41=dimm(mapda(0,3),sa3)
       nhelp42=dimm(mapda(0,4),sa4)
       if ((ntest2.eq.1).and.(sa3.eq.sa4)) then
       nhelp4=nhelp41*(nhelp41-1)/2
       else
       nhelp4=nhelp41*nhelp42
       end if
c
       mvec(ix,1)=nhelp1
       mvec(ix,2)=mapda(ia,1)
       mvec(ix,3)=mapdb(ib,1)
       mvec(ix,4)=mapdc(iy,1)
       mvec(ix,5)=nhelp2
       mvec(ix,6)=nhelp4
       mvec(ix,7)=0
c
       ix=ix+1
c
 50     continue
 100    continue
       ix=ix-1
c
       return
       end
c
c     --------------------------------
c
       subroutine cct3_grc32y (mapda,mapdb,mapdc,mapia,mapib,mapic,
     & mvec,ssa,ssb,possc0,ix)
c
#include "t31.fh"
c
       integer mapda(0:512,1:6)
       integer mapdb(0:512,1:6)
       integer mapdc(0:512,1:6)
c
       integer mapia(1:8,1:8,1:8)
       integer mapib(1:8,1:8,1:8)
       integer mapic(1:8,1:8,1:8)
c
       integer mvec(1:4096,1:7)
       integer possc0
       integer ssa,ssb
c
c     help variables
c
       integer nhelp1,nhelp2,nhelp4
       integer nhelp41,nhelp42
       integer ntest2
       integer sa1,sa2,sa3,sb1,sb2,sa23
       integer ia,ib,iy,ix
       integer possct
c
c     sctructure A(p,qr)*B(qr)=YC(p)
c
c1.0  prepare mapdc,mapic
c
       call cct3_grc0 (1,0,mapda(0,1),0,0,0,mmul(ssa,ssb),
     & possc0,possct,mapdc,mapic)
c
c1.1  define limitations - A p,q>r must be tested - ntest1
c     - B p>q must be tested - ntest2
c
c      if (mapda(0,6).eq.2) then
c      ntest1=1
c      else
c      ntest1=0
c      end if
c
       if (mapdb(0,6).eq.1) then
       ntest2=1
       else
       ntest2=0
       end if
c
c1.2  def symm states and test the limitations
c
       ix=1
       do 100 sb1=1,nsym
       sa2=sb1
c
       sb2=mmul(ssb,sb1)
       sa3=sb2
       sa23=mmul(sa2,sa3)
       if ((ntest2.eq.1).and.(sb1.lt.sb2)) then
c     Meggie out
       goto 100
       end if
c
       sa1=mmul(ssa,sa23)
c
c1.3  def mvec,mapdc and mapdi
c
       ia=mapia(sa1,sa2,1)
       ib=mapib(sb1,1,1)
       iy=mapic(1,1,1)
c
c     yes/no
       if ((mapda(ia,2).gt.0).and.(mapdb(ib,2).gt.0)) then
       nhelp1=1
       else
       goto 100
       end if
c
c     rowA
       nhelp2=dimm(mapda(0,1),sa1)
c
c     sum
       nhelp41=dimm(mapda(0,2),sa2)
       nhelp42=dimm(mapda(0,3),sa3)
       if ((ntest2.eq.1).and.(sa2.eq.sa3)) then
       nhelp4=nhelp41*(nhelp41-1)/2
       else
       nhelp4=nhelp41*nhelp42
       end if
c
       mvec(ix,1)=nhelp1
       mvec(ix,2)=mapda(ia,1)
       mvec(ix,3)=mapdb(ib,1)
       mvec(ix,4)=mapdc(iy,1)
       mvec(ix,5)=nhelp2
       mvec(ix,6)=nhelp4
       mvec(ix,7)=0
c
       ix=ix+1
c
 100    continue
       ix=ix-1
c
       return
       end
c
c     -----------------------------
c
       subroutine cct3_grc0 (nind,typ,typp,typq,typr,typs,stot,
     & poss0,posst,mapd,mapi)
c
c     this routine defines mapd and mapi for given intermediat
c
c     !N.B. (this routine cannot run with +OP2)
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
c
c     To fix some compiler warnings
      poss=0
      i=0
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
#ifdef _NOTUSED_
c
c     -----------------------------
c
       subroutine cct3_multdot (wrk,wrksize,
     & nind,mapda,mapia,ssa,mapdb,mapib,ssb,scalar,
     &                     rc)
c
c     this routine do dot product
c     scalar = A(ind)*B(ind)
c
c     nind        - number of indices in both A B (I)
c     mapda        - direct map of A (I)
c     mapia        - inverse map of A (I)
c     ssa        - symmetry state of A (I)
c     mapdb        - direct map of A (I)
c     mapib        - inverse map of A (I)
c     ssb        - symmetry state of B (I)
c     scalar        - final dot product (O)
c     rc        - return (error) code (O)
c
c     N.B. A and B must have same ordering of indexes
c
c     indA     indB     Implemented
c     1        1         Not yet
c     2        2           Yes
c     3        3           Yes
c     4        4           Yes
c     more              No
c
c
#include "wrk.fh"
c
       integer nind,ssa,ssb,rc
       real*8 scalar
c
       integer mapda(0:512,1:6)
       integer mapia(1:8,1:8,1:8)
c
       integer mapdb(0:512,1:6)
       integer mapib(1:8,1:8,1:8)
c
c     help variables
c
       integer symp,symq,symr,iia,iib,possa,possb,nhelp,length
       real*8 scal
c
       rc=0
c
cT    some tests
c
       do 10 nhelp=1,nind
       if (mapda(0,nhelp).ne.mapdb(0,nhelp)) then
c     RC =1 : nonidentical types of indices (NCI/Stup)
       rc=1
       return
       end if
 10     continue
c
       if (mapda(0,5).ne.mapdb(0,5)) then
c     RC =2 : nonidentical number of symmtry blocks in A and B (Stup)
       rc=2
       return
       end if
c
       if (mapda(0,6).ne.mapdb(0,6)) then
c     RC =3 : nonidentical type of A and B (Stup)
       rc=3
       return
       end if
c
       if (ssa.ne.ssb) then
c     RC =4 : nonidentical symmetry state of A and B (Stup)
       rc=4
       return
       end if
c

       if (nind.eq.4) then
c
cI    4 index matrices
c
       scalar=0.0d0
       do 100 iia=1,mapda(0,5)
c
cI.1  def parameters of A
       symp=mapda(iia,3)
       symq=mapda(iia,4)
       symr=mapda(iia,5)
c     syms is redundant
       possa=mapda(iia,1)
c
cI.2  def parameters of B
       iib=mapib(symp,symq,symr)
       possb=mapdb(iib,1)
c
cI.3  length must be common for both A and B
       length=mapda(iia,2)
c
       if (length.gt.0) then
       call cct3_mr0u3wt (length,length,length,1,1,wrk(possa),
     & wrk(possb),scal)
       scalar=scalar+scal
       end if
c
 100    continue
c
       else if (nind.eq.3) then
c
cII   3 index matrices
c
       scalar=0.0d0
       do 200 iia=1,mapda(0,5)
c
cII.1 def parameters of A
       symp=mapda(iia,3)
       symq=mapda(iia,4)
c     symr is redundant
       possa=mapda(iia,1)
c
cII.2 def parameters of B
       iib=mapib(symp,symq,1)
       possb=mapdb(iib,1)
c
cII.3 length must be common for both A and B
       length=mapda(iia,2)
c
       if (length.gt.0) then
       call cct3_mr0u3wt (length,length,length,1,1,wrk(possa),
     & wrk(possb),scal)
       scalar=scalar+scal
       end if
c
 200    continue
c
       else if (nind.eq.2) then
c
cIII  2 index matrices
c
       scalar=0.0d0
       do 300 iia=1,mapda(0,5)
c
cIII.1def parameters of A
       symp=mapda(iia,3)
c     symq is redundant
       possa=mapda(iia,1)
c
cIII.2def parameters of B
       iib=mapib(symp,1,1)
       possb=mapdb(iib,1)
c
cIII.3length must be common for both A and B
       length=mapda(iia,2)
c
       if (length.gt.0) then
       call cct3_mr0u3wt (length,length,length,1,1,wrk(possa),
     & wrk(possb),scal)
       scalar=scalar+scal
       end if
c
 300    continue
c
       else if (nind.eq.1) then
c
cIV   1 index matrices
c     RC=5 : 1 index matrices (NCI)
       rc=5
       return
c
       else
cV    more than 4 index matrices
c     RC=6 : more than 4 index matrices (NCI)
       rc=6
       return
       end if
c
c
       return
       end
c
c     -----------------------------
c
#endif
