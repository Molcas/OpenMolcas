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

subroutine ampack(wrk,wrksize,syma,symm,symp,symq,a,vint,ndimv1,ndimv2,ndimv3,ammap)
! this routine packs corresponding parts to ab direct acc. file
! from given integrals <_a,bb,p,q> read in vint
!
! syma  - irrep of first index (I)
! symm  - irrep of 2nd index - m (I)
! symp  - irrep of p (I)
! symq  - irrep of q (I)
! a     - pivot virtual index (counted in nvb set) (I)
! vint  - array of integrals <_a,mm,p,q> (I)
! ndimv1- 1st dimension of vint - norb(symb) (I)
! ndimv2- 2nd dimension of vint - norb(symp) (I)
! ndimv3- 3rd dimension of vint - norb(symq) (I)
! ammap - map for storing of addresses in DA file TEMPDA2 (I)

#include "wrk.fh"
#include "ccsort.fh"
#include "reorg.fh"
integer syma, symm, symp, symq, a, ndimv1, ndimv2, ndimv3
real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
integer ammap(1:mbas,1:8,1:8)
! help variables
integer m, p, q, pq, irec0, length

!T if there are no a, or no integrals in _a_mpq block return

if (nvb(syma)*noa(symm)*norb(symp)*norb(symq) == 0) then
  return
end if

! def length of _a(m,p,q) block

length = noa(symm)*norb(symp)*norb(symq)

! map _a(mpq) block into #v3

pq = poss30-1

do q=1,norb(symq)
  do p=1,norb(symp)
    do m=1,noa(symm)
      pq = pq+1
      wrk(pq) = vint(m,p,q)
    end do
  end do
end do

! put this block to appropriate position in direct access file

irec0 = ammap(a,symm,symp)
call dawrite(lunda2,irec0,wrk(poss30),length,recl)

return

end subroutine ampack
