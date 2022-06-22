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
       subroutine initintabc1
c     this routine write corresponding mapd and mapi to INTAB
c     for nonsymetrical (C1) case
c
#include "reorg.fh"
#include "ccsort.fh"
c
c     help variables
c
       integer nhelp,length,symp,symq,symab
       integer poss,ii,syma,symb,rc
c
c*    def symab
       syma=1
       symb=1
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
c
       return
       end
