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
       subroutine addpqij (wrk,wrksize,
     & symp,symq,symi,symj,p,vint,ndimv1,ndimv2,
     &                     ndimv3)
c
c     this routine add corresponding part to <pq|ij> record (#1)
c     comming from readed integrals with pivot index p vint_p(q,i,j)
c
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer symi,symj,symp,symq,p,ndimv1,ndimv2,ndimv3
       real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
c
c     help variables
c
       integer ii,ij,i,j,poss0,possij0,q,pqij
c
c*    find number of this symmetry combination
c     and initial possition of this symmetry block in (1)
c
       ii=mapi1(symp,symq,symi)
       poss0=mapd1(ii,1)
c
cT0   if symi<symj return
       if (symi.lt.symj) then
       return
       end if
c
cT1   return, if length is 0
       if (mapd1(ii,2).eq.0) then
       return
       end if
c
       do 1000 j=1,noa(symj)
       do 1001 i=1,noa(symi)
c
c*    def ij index and initial possition for <p,q,i,j> integral
c
       ij=(j-1)*noa(symi)+i
       possij0=poss0+(norb(symp)*norb(symq))*(ij-1)
c
       do 200 q=1,norb(symq)
       pqij=possij0-1+norb(symp)*(q-1)+p
       wrk(pqij)=vint(q,i,j)
 200    continue
c
 1001   continue
 1000   continue
c
       return
       end
