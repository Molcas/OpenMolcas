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
       subroutine abpack (wrk,wrksize,
     & syma,symb,symp,symq,a,vint,ndimv1,ndimv2,
     & ndimv3,abmap)
c
c     this routine pack corresponding parts to ab direct acc. file
c     from given integrals <_a,bb,p,q> readed in vint
c
c     syma  - irrep of first index (I)
c     symb  - irrep of 2.nd index (I)
c     symp  - irrep of p (I)
c     symq  - irrep of q (I)
c     a- pivot virtual index (counted in nvb set) (I)
c     vint  - array of integrals <_a,bb,p,q> (I)
c     ndimv1- 1.st dimension of vint - norb(symb) (I)
c     ndimv2- 2.nd dimension of vint - norb(symp) (I)
c     ndimv3- 3.rd dimension of vint - norb(symq) (I)
c     abmap - map for storing of addresses in DA file TEMPDA1 (I)
c
#include "wrk.fh"
#include "ccsort.fh"
#include "reorg.fh"
c
       integer syma,symb,symp,symq,a,ndimv1,ndimv2,ndimv3
       real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
       integer abmap(1:mbas,1:mbas,1:8)
c
c     help variables
c
       integer p,q,pq,irec0,length,b,bup,bvint
c
cT    if there are no ab pair, or no integrals in _a_bpq block return
c
       if (nvb(syma)*nvb(symb)*norb(symp)*norb(symq).eq.0) then
       return
       end if
c
c*    def length of _a_b(p,q) block
c
       length=norb(symp)*norb(symq)
c
       if (syma.eq.symb) then
       bup=a
       else
       bup=nvb(symb)
       end if
c
c*    cycle over b for given a
c
       do 1000 b=1,bup
       bvint=nob(symb)+b
c
c*    map _a_b(pq) block into #v3
c
       pq=poss30-1
       do 100 q=1,norb(symq)
       do 101 p=1,norb(symp)
       pq=pq+1
       wrk(pq)=vint(bvint,p,q)
 101    continue
 100    continue
c
c*    put this block to iappropriate possition in direct acces file
c
       irec0=abmap(a,b,symp)
       call dawrite (lunda1,irec0,wrk(poss30),length,recl)
c
 1000   continue
c
       return
       end
