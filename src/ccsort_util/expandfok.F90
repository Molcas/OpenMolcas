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
! it also defines new mapd2,mapi2

use ccsort_global, only: mapd2, mapi2, NORB, NSYM, pos20
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
real(kind=wp), intent(in) :: fok(*)
integer(kind=iwp) :: p, postemp, pqfok, pqwrk, q, qpwrk, symp

! set mapi zero

mapi2(1:nsym,1:nsym,1:nsym) = 0

! def zeroth row of mapd

mapd2(0,1) = 5
mapd2(0,2) = 5
mapd2(0,3) = 0
mapd2(0,4) = 0
mapd2(0,5) = nsym
mapd2(0,6) = 0

postemp = pos20
pqfok = 0
do symp=1,nsym

  ! def mapd,mapi

  mapd2(symp,1) = postemp
  mapd2(symp,2) = norb(symp)*norb(symp)
  mapd2(symp,3) = symp
  mapd2(symp,4) = symp
  mapd2(symp,5) = 1
  mapd2(symp,6) = 1
  mapi2(symp,1,1) = symp

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

  postemp = postemp+mapd2(symp,2)

end do

return

end subroutine expandfok
