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
       subroutine ampack (wrk,wrksize,
     & syma,symm,symp,symq,a,vint,ndimv1,ndimv2,ndimv3,
     & ammap)
c
c     this routine pack corresponding parts to ab direct acc. file
c     from given integrals <_a,bb,p,q> readed in vint
c
c     syma  - irrep of first index (I)
c     symm  - irrep of 2.nd index - m (I)
c     symp  - irrep of p (I)
c     symq  - irrep of q (I)
c     a- pivot virtual index (counted in nvb set) (I)
c     vint  - array of integrals <_a,mm,p,q> (I)
c     ndimv1- 1.st dimension of vint - norb(symb) (I)
c     ndimv2- 2.nd dimension of vint - norb(symp) (I)
c     ndimv3- 3.rd dimension of vint - norb(symq) (I)
c     ammap - map for storing of addresses in DA file TEMPDA2 (I)
c
#include "wrk.fh"
#include "ccsort.fh"
#include "reorg.fh"
c
       integer syma,symm,symp,symq,a,ndimv1,ndimv2,ndimv3
       real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
       integer ammap(1:mbas,1:8,1:8)
c
c     help variables
c
       integer m,p,q,pq,irec0,length
c
cT    if there are no a, or no integrals in _a_mpq block return
c
       if (nvb(syma)*noa(symm)*norb(symp)*norb(symq).eq.0) then
       return
       end if
c
c*    def length of _a(m,p,q) block
c
       length=noa(symm)*norb(symp)*norb(symq)
c
c*    map _a(mpq) block into #v3
c
       pq=poss30-1
c
       do 100 q=1,norb(symq)
       do 101 p=1,norb(symp)
       do 102 m=1,noa(symm)
       pq=pq+1
       wrk(pq)=vint(m,p,q)
 102    continue
 101    continue
 100    continue
c
c*    put this block to iappropriate possition in direct acces file
c
       irec0=ammap(a,symm,symp)
       call dawrite (lunda2,irec0,wrk(poss30),length,recl)
c
       return
       end
