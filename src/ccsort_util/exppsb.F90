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
       subroutine exppsb (symp,symq,symr,syms,                          &
     &                    valn,jn,kn,ln)
!
!     this routine realize expansion of symmetry block
!     symp,symq,symr,syms <p,q|r,s>  <-> (IJ|KL), provided such integrals exists
!     It found corresponding (IJ|KL) and expand it to opened
!     NORB(symp) TEMP files with a structure
!     indq,indr,inds,value, each TEMP for one p
!
!     N.B. This process can be accelerated, if exppbs would be
!     divided into exppsb1-8, each for given typ, since this
!     routine is common for all types.
!
!     types of (ij|kl) NI,J,K,L defined in III
!
!     1 - si=sk, si=sj, sk=sl
!     2 - si=sk, si=sj, sk>sl
!     3 - si=sk, si>sj, sk=sl
!     4 - si=sk, si>sj, sk>sl
!     5 - si>sk, si=sj, sk=sl
!     6 - si>sk, si=sj, sk>sl
!     7 - si>sk, si>sj, sk=sl
!     8 - si>sk, si>sj, sk>sl
!
!
       implicit real*8 (a-h,o-z)


#include "SysDef.fh"
#include "reorg.fh"
#include "ccsort.fh"
!
       integer symp,symq,symr,syms
       real*8 valn(1:nsize,1:mbas)
       integer jn(1:nsize,1:mbas)
       integer kn(1:nsize,1:mbas)
       integer ln(1:nsize,1:mbas)
!
!     help variables
!
       integer idis13,indtemp
       integer ni,nj,nk,nl,nsi,nsj,nsk,nsl,i1,j1,k1,l1
       integer iup,ilow,jup,jlow,kup,lup,iold,jold,kold,lold
!
       integer nhelp1,nhelp2,m3
       integer yes234,yes5,yes678
       integer typp
       integer ind(1:4)
#include "tratoc.fh"
       integer INDMAX
       parameter (INDMAX=nTraBuf)
       REAL*8    TWO(INDMAX)
!
!I    get  adress
       idis13=idis(symp,symq,symr)
!
!II   prepairing nshow vector
!
       do nhelp1=1,norb(symp)
       nshow(nhelp1)=0
       end do
!
!III.1define order of indices
!
       ni=np(symp,symq,symr)
       nj=nq(symp,symq,symr)
       nk=nr(symp,symq,symr)
       nl=ns(symp,symq,symr)
!
!III.2def yes1-8
!
       typp=typ(symp,symq,symr)
!
!:1   combination (ij|kl) -> (ij|kl)
!     used in types: 1,2,3,4,5,6,7,8 (all)
!      yes1=1
!
!:2   combination (ij|kl) -> (ji|kl)
!:3   combination (ij|kl) -> (ij|lk)
!:4   combination (ij|kl) -> (ji|lk)
!     used in types: 1,5 since 2,3,6,7 never appear
       if ((typp.eq.1).or.(typp.eq.5)) then
       yes234=1
       else
       yes234=0
       end if
!
!:5   combination (ij|kl) -> (kl|ij)
!     used in types: 1,2,3,4
       if ((typp.ge.1).and.(typp.le.4)) then
       yes5=1
       else
       yes5=0
       end if
!
!:6   combination (ij|kl) -> (lk|ij)
!:7   combination (ij|kl) -> (kl|ji)
!:8   combination (ij|kl) -> (lk|ji)
!     used in types: 1 (since 2,3 never appeard)
       if (typp.eq.1) then
       yes678=1
       else
       yes678=0
       end if
!
!
!     define NSI,NSJ,NSK,NSL
       ind(ni)=symp
       ind(nj)=symq
       ind(nk)=symr
       ind(nl)=syms
       NSI=ind(1)
       NSJ=ind(2)
       NSK=ind(3)
       NSL=ind(4)
!
       indtemp=indmax+1
       KUP=NORB(NSK)
       DO 401 KOLD=1,KUP
         if (fullprint.ge.3) write (6,*) ' * K ind ',KOLD
!
       LUP=NORB(NSL)
       IF (NSK.EQ.NSL) LUP=KOLD
       DO 402 LOLD=1,LUP
         if (fullprint.ge.3) write (6,*) ' ** L ind ',LOLD
!
       ILOW=1
       IF (NSI.EQ.NSK) ILOW=KOLD
       IUP=NORB(NSI)
       DO 403 IOLD=ILOW,IUP
         if (fullprint.ge.3) write (6,*) ' *** I ind ',IOLD
!
       JLOW=1
       IF (NSI.EQ.NSK.AND.IOLD.EQ.KOLD) JLOW=LOLD
       JUP=NORB(NSJ)
       IF (NSI.EQ.NSJ) JUP=IOLD
       DO 404 JOLD=JLOW,JUP
         if (fullprint.ge.3) write (6,*) ' **** J ind ',JOLD
!
!
!*    read block of integrals if necessary
!
       if (indtemp.eq.(indmax+1)) then
       indtemp=1
!     read block
       CALL dDAFILE(LUINTM,2,TWO,INDMAX,IDIS13)
       end if
!
!*    write integrals to appropriate positions
!
       val1=TWO(indtemp)
!
!:1   combination (ij|kl) -> (ij|kl)
!     since yes1 is always 1, if structure is skipped
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
!
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
!
       if (m3.eq.nsize) then
       call zasun (i1,nsize,                                            &
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
!
 21     if (yes234.eq.1) then
!
!:2   combination (ij|kl) -> (ji|kl)
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
!
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
!
       if (m3.eq.nsize) then
       call zasun (i1,nsize,                                            &
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
!
!:3   combination (ij|kl) -> (ij|lk)
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
!
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
!
       if (m3.eq.nsize) then
       call zasun (i1,nsize,                                            &
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
!
!:4   combination (ij|kl) -> (ji|lk)
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
!
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
!
       if (m3.eq.nsize) then
       call zasun (i1,nsize,                                            &
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
!
       end if
!
!:5   combination (ij|kl) -> (kl|ij)
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
!
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
!
       if (m3.eq.nsize) then
       call zasun (i1,nsize,                                            &
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
       end if
!
 61     if (yes678.eq.1) then
!
!:6   combination (ij|kl) -> (lk|ij)
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
!
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
!
       if (m3.eq.nsize) then
       call zasun (i1,nsize,                                            &
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
!
!:7   combination (ij|kl) -> (kl|ji)
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
!
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
!
       if (m3.eq.nsize) then
       call zasun (i1,nsize,                                            &
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
!
!:8   combination (ij|kl) -> (lk|ji)
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
!
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
!
       if (m3.eq.nsize) then
       call zasun (i1,nsize,                                            &
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
!
       end if
!
 100    indtemp=indtemp+1
!
 404    CONTINUE
 403    CONTINUE
 402    CONTINUE
 401    CONTINUE
!
!
!IV   write the rest integrals if needed
!
       do nhelp1=1,norb(symp)
       nhelp2=nshow(nhelp1)
       if (nhelp2.gt.0) then
       call zasun (nhelp1,nhelp2,                                       &
     &             valn,jn,kn,ln)
       end if
       end do
!
       return
       end
