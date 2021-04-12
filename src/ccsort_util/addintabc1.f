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
       subroutine addintabc1 (wrk,wrksize,
     & a,vint,ndimv)
c
c     this routine add integrals <_a,_b|p,q> for given a
c     for nonsymmetrical (C1) case
c     from integrals vv _a(u,p,q)
c
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer a,ndimv
       real*8 vint(1:ndimv,1:ndimv,1:ndimv)
c
c     help variables
c
       integer poss,b,bvint,p,q,length
c
c
cT    if there are no _a_b,pq integrals in this symab,
c     skip sumation over ab
c
       if (nvb(1).eq.0) then
       return
       end if
c
c*    loop over b
c
       do 1000 b=1,a
       bvint=b+nob(1)
c
c     map <_a,b|p,q> to wrk in #3
       poss=poss30
       do 1010 q=1,norb(1)
       do 1011 p=1,norb(1)
       wrk(poss)=vint(bvint,p,q)
       poss=poss+1
 1011   continue
 1010   continue
c
c
c**   since there must be some integrals, write them to TEMPAB
c
       length=poss-poss30
       call dawri (lunab,length,wrk(poss30))
c
 1000   continue
c
       return
       end
