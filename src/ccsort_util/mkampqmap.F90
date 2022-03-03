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

subroutine mkampqmap(ammap,syma,rc)
! this routine prepares ammap

use ccsort_global, only: noa, mbas, NORB, NSYM, nvb, reclen
use Symmetry_Info, only: Mul
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(_OUT_) :: ammap(mbas,8,8)
integer(kind=iwp), intent(in) :: syma
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: a, irec, lengthmpq, nrecc, nrest, symam, symm, symp, symq

rc = 0

!T test, if there are any a in this symmetry

if (nvb(syma) == 0) then
  rc = 1
  ! RC=1 : there are no a in this symmetry
  return
end if

! def initial address

irec = 1

! loop over all combinations

do symm=1,nsym
  symam = mul(syma,symm)
  do symp=1,nsym
    symq = mul(symam,symp)

    ! define number of records required to store this block
    ! and determine shift in initial positions

    lengthmpq = noa(symm)*norb(symp)*norb(symq)
    nrecc = int(lengthmpq/reclen)
    nrest = lengthmpq-nrecc*reclen
    if (nrest > 0) nrecc = nrecc+1

    do a=1,nvb(syma)

      ammap(a,symm,symp) = irec
      irec = irec+nrecc

    end do
  end do
end do

return

end subroutine mkampqmap
