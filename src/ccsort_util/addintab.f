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
       subroutine addintab (wrk,wrksize,
     & syma,symb,abmap)
c
c     this routine add contribution to opened INTAB1 file,
c     comming from ab syma,symb
c
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer syma,symb
       integer abmap(1:mbas,1:mbas,1:8)
c
c     help variables
c
       integer nhelp,length,symp,symq,symab,irec0,poss3
       integer poss,a,b,bup,ii,rc
c
c*    def symab
       symab=mul(syma,symb)
c
c*    make mapd3,mapi3 for <_a_b|pq>
c
c**   set mapi3=0 (partly)
c
       do 100 nhelp=1,nsym
       do 101 symq=1,nsym
       do 102 symp=1,nsym
       mapi3(symp,symq,nhelp)=0
 102    continue
 101    continue
 100    continue
c
c**   def 0-th row
c
       mapd3(0,1)=5
       mapd3(0,2)=5
       mapd3(0,3)=0
       mapd3(0,4)=0
       mapd3(0,5)=nsym
       mapd3(0,6)=0
c
c**   def other rows
c
       poss=poss30
       do 200 ii=1,nsym
c
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
c
 200    continue
c
c*    write mapd,mapi to INTAB
       call dawrtmap (lunab,mapd3,mapi3,rc)
c
cT    if there are no _a_b,pq integrals in this symab,
c     skip sumation over ab
c
       if ((mapd3(nsym,1)+mapd3(nsym,2)).eq.poss30) then
       return
       end if
c
c*    loop over a,b
c
       do 1000 a=1,nvb(syma)
c
       if (syma.eq.symb) then
       bup=a
       else
       bup=nvb(symb)
       end if
c
       do 1001 b=1,bup
c
c**   loop over symp
c
       do 500 symp=1,nsym
c
c***  def irec0 for this a,b,symp in TEMPDA1
       irec0=abmap(a,b,symp)
c
c***  def corresponding possition and length in #3
       ii=mapi3(symp,1,1)
       poss3=mapd3(ii,1)
       length=mapd3(ii,2)
c
c***  read this block to #3
       if (length.gt.0) then
       call daread (lunda1,irec0,wrk(poss3),length,recl)
       end if
c
 500    continue
c
c**   since there must be some integrals, write them to TEMPAB
c
       call deflength (mapd3,length)
       call dawri (lunab,length,wrk(poss30))
c
 1001   continue
 1000   continue
c
       return
       end
