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

subroutine addintabc1(wrk,wrksize,a,vint,ndimv)
! this routine adds integrals <_a,_b|p,q> for given a
! for nonsymmetrical (C1) case
! from integrals vv _a(u,p,q)

use ccsort_global, only: lunab, nob, NORB, nvb, pos30
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize, a, ndimv
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
real(kind=wp), intent(in) :: vint(ndimv,ndimv,ndimv)
integer(kind=iwp) :: b, bvint, length, p, pos, q

!T if there are no _a_b,pq integrals in this symab, skip summation over ab

if (nvb(1) == 0) return

! loop over b

do b=1,a
  bvint = b+nob(1)

  ! map <_a,b|p,q> to wrk in #3
  pos = pos30
  do q=1,norb(1)
    do p=1,norb(1)
      wrk(pos) = vint(bvint,p,q)
      pos = pos+1
    end do
  end do

  ! since there must be some integrals, write them to TEMPAB

  length = pos-pos30
  call dawri(lunab,length,wrk(pos30))

end do

return

end subroutine addintabc1
