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
       subroutine expandfok (wrk,wrksize,
     & fok)
c
c     This routine expand fok operator to #2
c     it also defines new mapd2,mapi2
c
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       real*8 fok(*)
c
c     help variables
c
       integer symp,symq,symr,posstemp,pqwrk,qpwrk,pqfok,p,q
c
c*    set mapi zero
c
       do 1 symr=1,nsym
       do 2 symq=1,nsym
       do 3 symp=1,nsym
       mapi2(symp,symq,symr)=0
 3      continue
 2      continue
 1      continue

c
c*    def zeroth row of mapd
c
       mapd2(0,1)=5
       mapd2(0,2)=5
       mapd2(0,3)=0
       mapd2(0,4)=0
       mapd2(0,5)=nsym
       mapd2(0,6)=0
c
       posstemp=poss20
       pqfok=0
       do 1000 symp=1,nsym
c
c*    def mapd,mapi
c
       mapd2(symp,1)=posstemp
       mapd2(symp,2)=norb(symp)*norb(symp)
       mapd2(symp,3)=symp
       mapd2(symp,4)=symp
       mapd2(symp,5)=1
       mapd2(symp,6)=1
       mapi2(symp,1,1)=symp
c
c*    expand
c
       do 100 p=1,norb(symp)
       do 101 q=1,p
c
c*    calc pq and qp possition in work and fok
c     and write integrals to this possitions
c
       pqwrk=posstemp+(norb(symp)*(p-1)+q)-1
       qpwrk=posstemp+(norb(symp)*(q-1)+p)-1
       pqfok=pqfok+1
       wrk(pqwrk)=fok(pqfok)
       wrk(qpwrk)=fok(pqfok)
c
 101    continue
 100    continue
c
       posstemp=posstemp+mapd2(symp,2)
c
 1000   continue
c
       return
       end
