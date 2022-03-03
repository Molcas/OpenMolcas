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

subroutine abpack(wrk,wrksize,syma,symb,symp,symq,a,vint,ndimv1,ndimv2,ndimv3,abmap)
! this routine packs corresponding parts to ab direct acc. file
! from given integrals <_a,bb,p,q> read in vint
!
! syma  - irrep of first index (I)
! symb  - irrep of 2nd index (I)
! symp  - irrep of p (I)
! symq  - irrep of q (I)
! a     - pivot virtual index (counted in nvb set) (I)
! vint  - array of integrals <_a,bb,p,q> (I)
! ndimv1- 1st dimension of vint - norb(symb) (I)
! ndimv2- 2nd dimension of vint - norb(symp) (I)
! ndimv3- 3rd dimension of vint - norb(symq) (I)
! abmap - map for storing of addresses in DA file TEMPDA1 (I)

use ccsort_global, only: lunda1, mbas, nob, NORB, nvb, pos30, reclen
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize, syma, symb, symp, symq, a, ndimv1, ndimv2, ndimv3, abmap(mbas,mbas,8)
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
real(kind=wp), intent(in) :: vint(ndimv1,ndimv2,ndimv3)
integer(kind=iwp) :: b, bup, bvint, irec0, length, p, pq, q

!T if there are no ab pair, or no integrals in _a_bpq block return

if (nvb(syma)*nvb(symb)*norb(symp)*norb(symq) == 0) return

! def length of _a_b(p,q) block

length = norb(symp)*norb(symq)

if (syma == symb) then
  bup = a
else
  bup = nvb(symb)
end if

! cycle over b for given a

do b=1,bup
  bvint = nob(symb)+b

  ! map _a_b(pq) block into #v3

  pq = pos30-1
  do q=1,norb(symq)
    do p=1,norb(symp)
      pq = pq+1
      wrk(pq) = vint(bvint,p,q)
    end do
  end do

  ! put this block to appropriate position in direct access file

  irec0 = abmap(a,b,symp)
  call dawrite(lunda1,irec0,wrk(pos30),length,reclen)

end do

return

end subroutine abpack
