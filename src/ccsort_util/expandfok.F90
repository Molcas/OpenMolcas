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

subroutine expandfok(wrk,wrksize,fok)
! This routine expands fok operator to #2
! it also defines new map2

use ccsort_global, only: map2, NORB, NSYM
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
real(kind=wp), intent(in) :: fok(*)
integer(kind=iwp) :: p, postemp, pqfok, pqwrk, q, qpwrk, symp

! set %i zero

map2%i(1:nsym,1:nsym,1:nsym) = 0

! def zeroth row of mapd

map2%d(0,1) = 5
map2%d(0,2) = 5
map2%d(0,3) = 0
map2%d(0,4) = 0
map2%d(0,5) = nsym
map2%d(0,6) = 0

postemp = map2%pos0
pqfok = 0
do symp=1,nsym

  ! def mapd,mapi

  map2%d(symp,1) = postemp
  map2%d(symp,2) = norb(symp)*norb(symp)
  map2%d(symp,3) = symp
  map2%d(symp,4) = symp
  map2%d(symp,5) = 1
  map2%d(symp,6) = 1
  map2%i(symp,1,1) = symp

  ! expand

  do p=1,norb(symp)
    do q=1,p

      ! calc pq and qp position in work and fok
      ! and write integrals to this positions

      pqwrk = postemp+(norb(symp)*(p-1)+q)-1
      qpwrk = postemp+(norb(symp)*(q-1)+p)-1
      pqfok = pqfok+1
      wrk(pqwrk) = fok(pqfok)
      wrk(qpwrk) = fok(pqfok)

    end do
  end do

  postemp = postemp+map2%d(symp,2)

end do

return

end subroutine expandfok
