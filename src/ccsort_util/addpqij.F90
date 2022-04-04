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

subroutine addpqij(wrk,wrksize,symp,symq,symi,symj,p,vint,ndimv1,ndimv2,ndimv3)
! this routine adds corresponding part to <pq|ij> record (#1)
! coming from read integrals with pivot index p vint_p(q,i,j)

use ccsort_global, only: mapd1, mapi1, noa, NORB
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize, symp, symq, symi, symj, p, ndimv1, ndimv2, ndimv3
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
real(kind=wp), intent(in) :: vint(ndimv1,ndimv2,ndimv3)
integer(kind=iwp) :: i, ii, ij, j, pos0, posij0, pqij, q

! find number of this symmetry combination
! and initial position of this symmetry block in (1)

ii = mapi1(symp,symq,symi)
pos0 = mapd1(ii,1)

!T0   if symi<symj return
if (symi < symj) return

!T1   return, if length is 0
if (mapd1(ii,2) == 0) return

do j=1,noa(symj)
  do i=1,noa(symi)

    ! def ij index and initial position for <p,q,i,j> integral

    ij = (j-1)*noa(symi)+i
    posij0 = pos0+(norb(symp)*norb(symq))*(ij-1)

    do q=1,norb(symq)
      pqij = posij0-1+norb(symp)*(q-1)+p
      wrk(pqij) = vint(q,i,j)
    end do

  end do
end do

return

end subroutine addpqij
