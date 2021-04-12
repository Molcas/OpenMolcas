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
       subroutine esb_ic_1 (symp,symq,symr,syms,
     c                      Vic,dimp,dimq,dimr,dims)
c
c     this routine realize expansion of symmetry block
c     symp,symq,symr,syms <p,q|r,s>  <-> (IJ|KL), provided such integrals exists
c     It found corresponding (IJ|KL) and expand it to
c     matrix vic (np,nq,nr,ns)
c

#include "SysDef.fh"
#include "reorg.fh"
#include "ccsort.fh"
c
        integer symp,symq,symr,syms
        integer dimp,dimq,dimr,dims
               real*8 Vic(1:dimp,1:dimq,1:dimr,1:dims)
        real*8 val1
c

       integer idis13,indtemp
       integer ni,nj,nk,nl,nsi,nsj,nsk,nsl,i1,j1,k1,l1
       integer iup,ilow,jup,jlow,kup,lup,iold,jold,kold,lold
c
c     help variables
       integer yes234,yes5,yes678
       integer typp
       integer ind(1:4)
#include "tratoc.fh"
       integer INDMAX
       parameter (INDMAX=nTraBuf)
       REAL*8    TWO(INDMAX)
c
cI    get  adress
       idis13=idis(symp,symq,symr)
c
cIII.1define order of indices
c
       ni=np(symp,symq,symr)
       nj=nq(symp,symq,symr)
       nk=nr(symp,symq,symr)
       nl=ns(symp,symq,symr)
c
cIII.2def yes1-8
c
       typp=typ(symp,symq,symr)
c
c:1   combination (ij|kl) -> (ij|kl)
c     used in types: 1,2,3,4,5,6,7,8 (all)
c      yes1=1
c
c:2   combination (ij|kl) -> (ji|kl)
c:3   combination (ij|kl) -> (ij|lk)
c:4   combination (ij|kl) -> (ji|lk)
c     used in types: 1,5 since 2,3,6,7 never appear
       if ((typp.eq.1).or.(typp.eq.5)) then
       yes234=1
       else
       yes234=0
       end if
c
c:5   combination (ij|kl) -> (kl|ij)
c     used in types: 1,2,3,4
       if ((typp.ge.1).and.(typp.le.4)) then
       yes5=1
       else
       yes5=0
       end if
c
c:6   combination (ij|kl) -> (lk|ij)
c:7   combination (ij|kl) -> (kl|ji)
c:8   combination (ij|kl) -> (lk|ji)
c     used in types: 1 (since 2,3 never appeard)
       if (typp.eq.1) then
       yes678=1
       else
       yes678=0
       end if
c
c
c     define NSI,NSJ,NSK,NSL
       ind(ni)=symp
       ind(nj)=symq
       ind(nk)=symr
       ind(nl)=syms
       NSI=ind(1)
       NSJ=ind(2)
       NSK=ind(3)
       NSL=ind(4)
C
       indtemp=indmax+1
       KUP=NORB(NSK)
       DO 401 KOLD=1,KUP
C
       LUP=NORB(NSL)
       IF (NSK.EQ.NSL) LUP=KOLD
       DO 402 LOLD=1,LUP
C
       ILOW=1
       IF (NSI.EQ.NSK) ILOW=KOLD
       IUP=NORB(NSI)
       DO 403 IOLD=ILOW,IUP
C
       JLOW=1
       IF (NSI.EQ.NSK.AND.IOLD.EQ.KOLD) JLOW=LOLD
       JUP=NORB(NSJ)
       IF (NSI.EQ.NSJ) JUP=IOLD
       DO 404 JOLD=JLOW,JUP
C
c
c*    read block of integrals if neccesarry
c
       if (indtemp.eq.(indmax+1)) then
       indtemp=1
c     read block
       CALL dDAFILE(LUINTM,2,TWO,INDMAX,IDIS13)
       end if
c
c*    write integrals to appropriate possitions
c
       val1=TWO(indtemp)
c
c:1   combination (ij|kl) -> (ij|kl)
c     since yes1 is always 1, if structure is skipped
       ind(1)=iold
       ind(2)=jold
       ind(3)=kold
       ind(4)=lold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimr) then
             if (l1.le.dims) then
             Vic(i1,j1,k1,l1)=val1
             end if
           end if
         end if
       end if
c
       if (yes234.eq.1) then
c
c:2   combination (ij|kl) -> (ji|kl)
       ind(1)=jold
       ind(2)=iold
       ind(3)=kold
       ind(4)=lold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimr) then
             if (l1.le.dims) then
             Vic(i1,j1,k1,l1)=val1
             end if
           end if
         end if
       end if
c
c:3   combination (ij|kl) -> (ij|lk)
       ind(1)=iold
       ind(2)=jold
       ind(3)=lold
       ind(4)=kold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimr) then
             if (l1.le.dims) then
             Vic(i1,j1,k1,l1)=val1
             end if
           end if
         end if
       end if
c
c:4   combination (ij|kl) -> (ji|lk)
       ind(1)=jold
       ind(2)=iold
       ind(3)=lold
       ind(4)=kold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimr) then
             if (l1.le.dims) then
             Vic(i1,j1,k1,l1)=val1
             end if
           end if
         end if
       end if
c
       end if
c
c:5   combination (ij|kl) -> (kl|ij)
       if (yes5.eq.1) then
       ind(1)=kold
       ind(2)=lold
       ind(3)=iold
       ind(4)=jold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimr) then
             if (l1.le.dims) then
             Vic(i1,j1,k1,l1)=val1
             end if
           end if
         end if
       end if
c
       end if
c
       if (yes678.eq.1) then
c
c:6   combination (ij|kl) -> (lk|ij)
       ind(1)=lold
       ind(2)=kold
       ind(3)=iold
       ind(4)=jold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimr) then
             if (l1.le.dims) then
             Vic(i1,j1,k1,l1)=val1
             end if
           end if
         end if
       end if
c
c:7   combination (ij|kl) -> (kl|ji)
       ind(1)=kold
       ind(2)=lold
       ind(3)=jold
       ind(4)=iold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimr) then
             if (l1.le.dims) then
             Vic(i1,j1,k1,l1)=val1
             end if
           end if
         end if
       end if
c
c:8   combination (ij|kl) -> (lk|ji)
       ind(1)=lold
       ind(2)=kold
       ind(3)=jold
       ind(4)=iold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimr) then
             if (l1.le.dims) then
             Vic(i1,j1,k1,l1)=val1
             end if
           end if
         end if
       end if
c
       end if
c
        indtemp=indtemp+1
c
 404    CONTINUE
 403    CONTINUE
 402    CONTINUE
 401    CONTINUE
C
c
       return
       end
