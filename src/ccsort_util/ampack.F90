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

use ccsort_global, only: lunda2, mbas, noa, NORB, nvb, pos30, reclen
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize, syma, symm, symp, symq, a, ndimv1, ndimv2, ndimv3, ammap(mbas,8,8)
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
real(kind=wp), intent(in) :: vint(ndimv1,ndimv2,ndimv3)
integer(kind=iwp) :: irec0, length, m, p, pq, q

!T if there are no a, or no integrals in _a_mpq block return

if (nvb(syma)*noa(symm)*norb(symp)*norb(symq) == 0) return

! def length of _a(m,p,q) block

length = noa(symm)*norb(symp)*norb(symq)

! map _a(mpq) block into #v3

pq = pos30-1

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
call dawrite(lunda2,irec0,wrk(pos30),length,reclen)

return

end subroutine ampack
