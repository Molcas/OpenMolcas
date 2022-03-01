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
       subroutine abpack (wrk,wrksize,                                  &
     & syma,symb,symp,symq,a,vint,ndimv1,ndimv2,                        &
     & ndimv3,abmap)
!
!     this routine pack corresponding parts to ab direct acc. file
!     from given integrals <_a,bb,p,q> readed in vint
!
!     syma  - irrep of first index (I)
!     symb  - irrep of 2.nd index (I)
!     symp  - irrep of p (I)
!     symq  - irrep of q (I)
!     a- pivot virtual index (counted in nvb set) (I)
!     vint  - array of integrals <_a,bb,p,q> (I)
!     ndimv1- 1.st dimension of vint - norb(symb) (I)
!     ndimv2- 2.nd dimension of vint - norb(symp) (I)
!     ndimv3- 3.rd dimension of vint - norb(symq) (I)
!     abmap - map for storing of addresses in DA file TEMPDA1 (I)
!
#include "wrk.fh"
#include "ccsort.fh"
#include "reorg.fh"
!
       integer syma,symb,symp,symq,a,ndimv1,ndimv2,ndimv3
       real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
       integer abmap(1:mbas,1:mbas,1:8)
!
!     help variables
!
       integer p,q,pq,irec0,length,b,bup,bvint
!
!T    if there are no ab pair, or no integrals in _a_bpq block return
!
       if (nvb(syma)*nvb(symb)*norb(symp)*norb(symq).eq.0) then
       return
       end if
!
!*    def length of _a_b(p,q) block
!
       length=norb(symp)*norb(symq)
!
       if (syma.eq.symb) then
       bup=a
       else
       bup=nvb(symb)
       end if
!
!*    cycle over b for given a
!
       do 1000 b=1,bup
       bvint=nob(symb)+b
!
!*    map _a_b(pq) block into #v3
!
       pq=poss30-1
       do 100 q=1,norb(symq)
       do 101 p=1,norb(symp)
       pq=pq+1
       wrk(pq)=vint(bvint,p,q)
 101    continue
 100    continue
!
!*    put this block to iappropriate possition in direct acces file
!
       irec0=abmap(a,b,symp)
       call dawrite (lunda1,irec0,wrk(poss30),length,recl)
!
 1000   continue
!
       return
       end
