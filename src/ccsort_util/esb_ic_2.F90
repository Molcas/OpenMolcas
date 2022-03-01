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
       subroutine esb_ic_2 (symp,symq,Vic,dimp,dimq,pqind)
!
!     this routine realize expansion of symmetry block
!     symp,symq,symr,syms <p,q|r,s>  <-> (IJ|KL),
!     (for case symp=symr,symq=syms)
!     provided such integrals exists
!     It found corresponding (IJ|KL) and expand it to
!     matrix vic (pr,qs)
!

#include "SysDef.fh"
#include "reorg.fh"
#include "ccsort.fh"
!
        integer symp,symq
        integer dimp,dimq
               real*8 Vic(1:(dimp*(dimp+1)/2),1:(dimq*(dimq+1)/2))
        real*8 val1
!

       integer idis13,indtemp
       integer ni,nj,nk,nl,nsi,nsj,nsk,nsl,i1,j1,k1,l1
       integer iup,ilow,jup,jlow,kup,lup,iold,jold,kold,lold
       integer pqind(1:mbas,1:mbas)
!
!     help variables
       integer i,j,maxx
       integer yes234,yes5,yes678
       integer typp
       integer ind(1:4)
#include "tratoc.fh"
       integer INDMAX
       parameter (INDMAX=nTraBuf)
       REAL*8    TWO(INDMAX)
!
!I        calc pqind
!
        if (dimp.ge.dimq) then
          maxx=dimp
        else
          maxx=dimq
        end if
!
        do i=1,maxx
        do j=1,maxx
          if (i.ge.j) then
            pqind(i,j)=i*(i-1)/2+j
          else
            pqind(i,j)=j*(j-1)/2+i
          end if
        end do
        end do
!
!II    get  adress
       idis13=idis(symp,symq,symp)
!
!III.1define order of indices
!
       ni=np(symp,symq,symp)
       nj=nq(symp,symq,symp)
       nk=nr(symp,symq,symp)
       nl=ns(symp,symq,symp)
!
!III.2def yes1-8
!
       typp=typ(symp,symq,symp)
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
       ind(nk)=symp
       ind(nl)=symq
       NSI=ind(1)
       NSJ=ind(2)
       NSK=ind(3)
       NSL=ind(4)
!
       indtemp=indmax+1
       KUP=NORB(NSK)
       DO 401 KOLD=1,KUP
!
       LUP=NORB(NSL)
       IF (NSK.EQ.NSL) LUP=KOLD
       DO 402 LOLD=1,LUP
!
       ILOW=1
       IF (NSI.EQ.NSK) ILOW=KOLD
       IUP=NORB(NSI)
       DO 403 IOLD=ILOW,IUP
!
       JLOW=1
       IF (NSI.EQ.NSK.AND.IOLD.EQ.KOLD) JLOW=LOLD
       JUP=NORB(NSJ)
       IF (NSI.EQ.NSJ) JUP=IOLD
       DO 404 JOLD=JLOW,JUP
!
!
!*    read block of integrals if neccesarry
!
       if (indtemp.eq.(indmax+1)) then
       indtemp=1
!     read block
       CALL dDAFILE(LUINTM,2,TWO,INDMAX,IDIS13)
       end if
!
!*    write integrals to appropriate possitions
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
       i1=ind(ni)
       k1=ind(nk)
!
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimp) then
             if (l1.le.dimq) then
             Vic(pqind(i1,k1),pqind(j1,l1))=val1
             end if
           end if
         end if
       end if
!
       if (yes234.eq.1) then
!
!:2   combination (ij|kl) -> (ji|kl)
       ind(1)=jold
       ind(2)=iold
       ind(3)=kold
       ind(4)=lold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
!
!
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimp) then
             if (l1.le.dimq) then
             Vic(pqind(i1,k1),pqind(j1,l1))=val1
             end if
           end if
         end if
       end if
!
!:3   combination (ij|kl) -> (ij|lk)
       ind(1)=iold
       ind(2)=jold
       ind(3)=lold
       ind(4)=kold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
!
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimp) then
             if (l1.le.dimq) then
             Vic(pqind(i1,k1),pqind(j1,l1))=val1
             end if
           end if
         end if
       end if
!
!:4   combination (ij|kl) -> (ji|lk)
       ind(1)=jold
       ind(2)=iold
       ind(3)=lold
       ind(4)=kold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
!
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimp) then
             if (l1.le.dimq) then
             Vic(pqind(i1,k1),pqind(j1,l1))=val1
             end if
           end if
         end if
       end if
!
       end if
!
!:5   combination (ij|kl) -> (kl|ij)
       if (yes5.eq.1) then
       ind(1)=kold
       ind(2)=lold
       ind(3)=iold
       ind(4)=jold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
!
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimp) then
             if (l1.le.dimq) then
             Vic(pqind(i1,k1),pqind(j1,l1))=val1
             end if
           end if
         end if
       end if
!
       end if
!
       if (yes678.eq.1) then
!
!:6   combination (ij|kl) -> (lk|ij)
       ind(1)=lold
       ind(2)=kold
       ind(3)=iold
       ind(4)=jold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
!
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimp) then
             if (l1.le.dimq) then
             Vic(pqind(i1,k1),pqind(j1,l1))=val1
             end if
           end if
         end if
       end if
!
!:7   combination (ij|kl) -> (kl|ji)
       ind(1)=kold
       ind(2)=lold
       ind(3)=jold
       ind(4)=iold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
!
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimp) then
             if (l1.le.dimq) then
             Vic(pqind(i1,k1),pqind(j1,l1))=val1
             end if
           end if
         end if
       end if
!
!:8   combination (ij|kl) -> (lk|ji)
       ind(1)=lold
       ind(2)=kold
       ind(3)=jold
       ind(4)=iold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
!
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimp) then
             if (l1.le.dimq) then
             Vic(pqind(i1,k1),pqind(j1,l1))=val1
             end if
           end if
         end if
       end if
!
       end if
!
        indtemp=indtemp+1
!
 404    CONTINUE
 403    CONTINUE
 402    CONTINUE
 401    CONTINUE
!
!
       return
       end
