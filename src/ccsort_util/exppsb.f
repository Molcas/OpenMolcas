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
       subroutine exppsb (symp,symq,symr,syms,
     &                    valn,jn,kn,ln)
c
c     this routine realize expansion of symmetry block
c     symp,symq,symr,syms <p,q|r,s>  <-> (IJ|KL), provided such integrals exists
c     It found corresponding (IJ|KL) and expand it to opened
c     NORB(symp) TEMP files with a structure
c     indq,indr,inds,value, each TEMP for one p
c
c     N.B. This process can be accelerated, if exppbs would be
c     divided into exppsb1-8, each for given typ, since this
c     routine is common for all types.
c
c     types of (ij|kl) NI,J,K,L defined in III
c
c     1 - si=sk, si=sj, sk=sl
c     2 - si=sk, si=sj, sk>sl
c     3 - si=sk, si>sj, sk=sl
c     4 - si=sk, si>sj, sk>sl
c     5 - si>sk, si=sj, sk=sl
c     6 - si>sk, si=sj, sk>sl
c     7 - si>sk, si>sj, sk=sl
c     8 - si>sk, si>sj, sk>sl
c
c
       implicit real*8 (a-h,o-z)


#include "SysDef.fh"
#include "reorg.fh"
#include "ccsort.fh"
c
       integer symp,symq,symr,syms
       real*8 valn(1:nsize,1:mbas)
       integer jn(1:nsize,1:mbas)
       integer kn(1:nsize,1:mbas)
       integer ln(1:nsize,1:mbas)
c
c     help variables
c
       integer idis13,indtemp
       integer ni,nj,nk,nl,nsi,nsj,nsk,nsl,i1,j1,k1,l1
       integer iup,ilow,jup,jlow,kup,lup,iold,jold,kold,lold
c
       integer nhelp1,nhelp2,m3
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
cII   prepairing nshow vector
c
       do nhelp1=1,norb(symp)
       nshow(nhelp1)=0
       end do
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
         if (fullprint.ge.3) write (6,*) ' * K ind ',KOLD
C
       LUP=NORB(NSL)
       IF (NSK.EQ.NSL) LUP=KOLD
       DO 402 LOLD=1,LUP
         if (fullprint.ge.3) write (6,*) ' ** L ind ',LOLD
C
       ILOW=1
       IF (NSI.EQ.NSK) ILOW=KOLD
       IUP=NORB(NSI)
       DO 403 IOLD=ILOW,IUP
         if (fullprint.ge.3) write (6,*) ' *** I ind ',IOLD
C
       JLOW=1
       IF (NSI.EQ.NSK.AND.IOLD.EQ.KOLD) JLOW=LOLD
       JUP=NORB(NSJ)
       IF (NSI.EQ.NSJ) JUP=IOLD
       DO 404 JOLD=JLOW,JUP
         if (fullprint.ge.3) write (6,*) ' **** J ind ',JOLD
C
c
c*    read block of integrals if necessary
c
       if (indtemp.eq.(indmax+1)) then
       indtemp=1
c     read block
       CALL dDAFILE(LUINTM,2,TWO,INDMAX,IDIS13)
       end if
c
c*    write integrals to appropriate positions
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
       if (symq.eq.syms) then
       if (l1.gt.j1) then
       goto 21
       end if
       end if
       i1=ind(ni)
       k1=ind(nk)
c
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
c
       if (m3.eq.nsize) then
       call zasun (i1,nsize,
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
c
 21     if (yes234.eq.1) then
c
c:2   combination (ij|kl) -> (ji|kl)
       ind(1)=jold
       ind(2)=iold
       ind(3)=kold
       ind(4)=lold
       j1=ind(nj)
       l1=ind(nl)
       if (symq.eq.syms) then
       if (l1.gt.j1) then
       goto 31
       end if
       end if
       i1=ind(ni)
       k1=ind(nk)
c
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
c
       if (m3.eq.nsize) then
       call zasun (i1,nsize,
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
c
c:3   combination (ij|kl) -> (ij|lk)
 31     ind(1)=iold
       ind(2)=jold
       ind(3)=lold
       ind(4)=kold
       j1=ind(nj)
       l1=ind(nl)
       if (symq.eq.syms) then
       if (l1.gt.j1) then
       goto 41
       end if
       end if
       i1=ind(ni)
       k1=ind(nk)
c
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
c
       if (m3.eq.nsize) then
       call zasun (i1,nsize,
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
c
c:4   combination (ij|kl) -> (ji|lk)
 41     ind(1)=jold
       ind(2)=iold
       ind(3)=lold
       ind(4)=kold
       j1=ind(nj)
       l1=ind(nl)
       if (symq.eq.syms) then
       if (l1.gt.j1) then
       goto 51
       end if
       end if
       i1=ind(ni)
       k1=ind(nk)
c
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
c
       if (m3.eq.nsize) then
       call zasun (i1,nsize,
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
c
       end if
c
c:5   combination (ij|kl) -> (kl|ij)
 51     if (yes5.eq.1) then
       ind(1)=kold
       ind(2)=lold
       ind(3)=iold
       ind(4)=jold
       j1=ind(nj)
       l1=ind(nl)
       if (symq.eq.syms) then
       if (l1.gt.j1) then
       goto 61
       end if
       end if
       i1=ind(ni)
       k1=ind(nk)
c
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
c
       if (m3.eq.nsize) then
       call zasun (i1,nsize,
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
       end if
c
 61     if (yes678.eq.1) then
c
c:6   combination (ij|kl) -> (lk|ij)
       ind(1)=lold
       ind(2)=kold
       ind(3)=iold
       ind(4)=jold
       j1=ind(nj)
       l1=ind(nl)
       if (symq.eq.syms) then
       if (l1.gt.j1) then
       goto 71
       end if
       end if
       i1=ind(ni)
       k1=ind(nk)
c
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
c
       if (m3.eq.nsize) then
       call zasun (i1,nsize,
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
c
c:7   combination (ij|kl) -> (kl|ji)
 71     ind(1)=kold
       ind(2)=lold
       ind(3)=jold
       ind(4)=iold
       j1=ind(nj)
       l1=ind(nl)
       if (symq.eq.syms) then
       if (l1.gt.j1) then
       goto 81
       end if
       end if
       i1=ind(ni)
       k1=ind(nk)
c
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
c
       if (m3.eq.nsize) then
       call zasun (i1,nsize,
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
c
c:8   combination (ij|kl) -> (lk|ji)
 81     ind(1)=lold
       ind(2)=kold
       ind(3)=jold
       ind(4)=iold
       j1=ind(nj)
       l1=ind(nl)
       if (symq.eq.syms) then
       if (l1.gt.j1) then
       goto 100
       end if
       end if
       i1=ind(ni)
       k1=ind(nk)
c
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
c
       if (m3.eq.nsize) then
       call zasun (i1,nsize,
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
c
       end if
c
 100    indtemp=indtemp+1
c
 404    CONTINUE
 403    CONTINUE
 402    CONTINUE
 401    CONTINUE
C
c
cIV   write the rest integrals if needed
c
       do nhelp1=1,norb(symp)
       nhelp2=nshow(nhelp1)
       if (nhelp2.gt.0) then
       call zasun (nhelp1,nhelp2,
     &             valn,jn,kn,ln)
       end if
       end do
c
       return
       end
