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
       subroutine addintab (wrk,wrksize,                                &
     & syma,symb,abmap)
!
!     this routine add contribution to opened INTAB1 file,
!     comming from ab syma,symb
!
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer syma,symb
       integer abmap(1:mbas,1:mbas,1:8)
!
!     help variables
!
       integer nhelp,length,symp,symq,symab,irec0,poss3
       integer poss,a,b,bup,ii,rc
!
!*    def symab
       symab=mul(syma,symb)
!
!*    make mapd3,mapi3 for <_a_b|pq>
!
!**   set mapi3=0 (partly)
!
       do 100 nhelp=1,nsym
       do 101 symq=1,nsym
       do 102 symp=1,nsym
       mapi3(symp,symq,nhelp)=0
 102    continue
 101    continue
 100    continue
!
!**   def 0-th row
!
       mapd3(0,1)=5
       mapd3(0,2)=5
       mapd3(0,3)=0
       mapd3(0,4)=0
       mapd3(0,5)=nsym
       mapd3(0,6)=0
!
!**   def other rows
!
       poss=poss30
       do 200 ii=1,nsym
!
       symp=ii
       symq=mul(symab,symp)
       length=norb(symp)*norb(symq)
       mapd3(ii,1)=poss
       mapd3(ii,2)=length
       mapd3(ii,3)=symp
       mapd3(ii,4)=symq
       mapd3(ii,5)=1
       mapd3(ii,6)=1
       mapi3(symp,1,1)=ii
       poss=poss+length
!
 200    continue
!
!*    write mapd,mapi to INTAB
       call dawrtmap (lunab,mapd3,mapi3,rc)
!
!T    if there are no _a_b,pq integrals in this symab,
!     skip sumation over ab
!
       if ((mapd3(nsym,1)+mapd3(nsym,2)).eq.poss30) then
       return
       end if
!
!*    loop over a,b
!
       do 1000 a=1,nvb(syma)
!
       if (syma.eq.symb) then
       bup=a
       else
       bup=nvb(symb)
       end if
!
       do 1001 b=1,bup
!
!**   loop over symp
!
       do 500 symp=1,nsym
!
!***  def irec0 for this a,b,symp in TEMPDA1
       irec0=abmap(a,b,symp)
!
!***  def corresponding possition and length in #3
       ii=mapi3(symp,1,1)
       poss3=mapd3(ii,1)
       length=mapd3(ii,2)
!
!***  read this block to #3
       if (length.gt.0) then
       call daread (lunda1,irec0,wrk(poss3),length,recl)
       end if
!
 500    continue
!
!**   since there must be some integrals, write them to TEMPAB
!
       call deflength (mapd3,length)
       call dawri (lunab,length,wrk(poss30))
!
 1001   continue
 1000   continue
!
       return
       end
