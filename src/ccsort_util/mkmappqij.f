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
       subroutine mkmappqij
c
c     this routine prepair mapd,mapi
c     for <pq|ij> for p,q, i>=j to mapd1,mapi1
c
#include "ccsort.fh"
#include "reorg.fh"
c
c     help variables
c
       integer symi,symj,symp,symq,sympq,sympqi
       integer nhelp,possition,length
c
c*    set mapi1 to zero
c
       do 1 symi=1,nsym
       do 2 symq=1,nsym
       do 3 symp=1,nsym
       mapi1(symp,symq,symi)=0
 3      continue
 2      continue
 1      continue
c
c     def zero-th row
c
       mapd1(0,1)=5
       mapd1(0,2)=5
       mapd1(0,3)=1
       mapd1(0,4)=1
       mapd1(0,6)=3

       nhelp=0
       possition=poss10
       do 100 symp=1,nsym
       do 101 symq=1,nsym
       sympq=mul(symp,symq)
       do 102 symi=1,nsym
       sympqi=mul(sympq,symi)
       symj=sympqi
       if (symj.gt.symi) goto 102
       nhelp=nhelp+1
c
c     calc. length
       length=noa(symi)*noa(symj)*NORB(symp)*NORB(symq)
c
       mapd1(nhelp,1)=possition
       mapd1(nhelp,2)=length
       mapd1(nhelp,3)=symp
       mapd1(nhelp,4)=symq
       mapd1(nhelp,5)=symi
       mapd1(nhelp,6)=symj
       possition=possition+length
c
       mapi1(symp,symq,symi)=nhelp
c
 102    continue
 101    continue
 100    continue
c
       mapd1(0,5)=nhelp
c
       return
       end
