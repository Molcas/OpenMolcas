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
       subroutine mkmapampq (syma)
c
c     this routine prepair mapd,mapi
c     for <am|pq> for given syma, m, p,q to mapd2,mapi2
c
#include "ccsort.fh"
#include "reorg.fh"
       integer syma
c
c     help variables
c
       integer symm,symp,symq,symmp
       integer nhelp,possition,length
c
c*    set mapi1 to zero
c
       do 1 symq=1,nsym
       do 2 symp=1,nsym
       do 3 symm=1,nsym
       mapi2(symm,symp,symq)=0
 3      continue
 2      continue
 1      continue
c
c     def zero-th row
c
       mapd2(0,1)=1
       mapd2(0,2)=5
       mapd2(0,3)=5
       mapd2(0,4)=0
       mapd2(0,6)=0

       nhelp=0
       possition=poss20
       do 100 symm=1,nsym
       do 101 symp=1,nsym
       symmp=mul(symm,symp)
       symq=mul(syma,symmp)
       nhelp=nhelp+1
c
c     calc. length
       length=noa(symm)*NORB(symp)*NORB(symq)
c
       mapd2(nhelp,1)=possition
       mapd2(nhelp,2)=length
       mapd2(nhelp,3)=symm
       mapd2(nhelp,4)=symp
       mapd2(nhelp,5)=symq
       mapd2(nhelp,6)=1
       possition=possition+length
c
       mapi2(symm,symp,1)=nhelp
c
 101    continue
 100    continue
c
       mapd2(0,5)=nhelp
c
       return
       end
