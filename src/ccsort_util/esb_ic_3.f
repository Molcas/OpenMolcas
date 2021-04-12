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
       subroutine esb_ic_3 (symp,Vic,dimp,pqind)
c
c     this routine realize expansion of symmetry block
c     symp,symq,symr,syms <p,q|r,s>  <-> (IJ|KL),
c     (for case symp=symr=symq=syms)
c     provided such integrals exists
c     It found corresponding (IJ|KL) and expand it to
c     matrix vic (prqs)
c

#include "SysDef.fh"
#include "reorg.fh"
#include "ccsort.fh"
c
           integer symp,dimp
CLD           integer symp,dimp,dimq
               real*8 Vic(1:(dimp*(dimp+1)/2)*((dimp*(dimp+1)/2)+1)/2)
        real*8 val1
c

       integer idis13,indtemp
       integer ni,nj,nk,nl,nsi,nsj,nsk,nsl,i1,j1,k1,l1
       integer iup,ilow,jup,jlow,kup,lup,iold,jold,kold,lold
       integer pqind(1:mbas,1:mbas)
c
c     help variables
       integer i,j,maxx,ik,jl,ikjl
       integer ind(1:4)
#include "tratoc.fh"
       integer INDMAX
       parameter (INDMAX=nTraBuf)
       REAL*8    TWO(INDMAX)
c
cI        calc pqind
c
CLD        if (dimp.ge.dimp) then
          maxx=dimp
CLD        else
CLD          maxx=dimq
CLD        end if
c
        do i=1,maxx
        do j=1,maxx
          if (i.ge.j) then
            pqind(i,j)=i*(i-1)/2+j
          else
            pqind(i,j)=j*(j-1)/2+i
          end if
        end do
        end do
c
cII    get  adress
       idis13=idis(symp,symp,symp)
c
cIII.1define order of indices
c
       ni=np(symp,symp,symp)
       nj=nq(symp,symp,symp)
       nk=nr(symp,symp,symp)
       nl=ns(symp,symp,symp)
c
c
c     define NSI,NSJ,NSK,NSL
       ind(ni)=symp
       ind(nj)=symp
       ind(nk)=symp
       ind(nl)=symp
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
c:2    def iklj
       ik=pqind(i1,k1)
       jl=pqind(j1,l1)
       if (ik.ge.jl) then
         ikjl=ik*(ik-1)/2+jl
       else
         ikjl=jl*(jl-1)/2+ik
       end if
c
       Vic(ikjl)=val1
c
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
