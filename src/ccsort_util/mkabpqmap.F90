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

subroutine mkabpqmap(abmap,syma,symb,rc)
! this routine prepares abmap

use ccsort_global, only: mbas, NORB, NSYM, nvb, reclen
use Symmetry_Info, only: Mul
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(_OUT_) :: abmap(mbas,mbas,8)
integer(kind=iwp), intent(in) :: syma, symb
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: a, b, bup, irec, lengthpq, nrecc, nrest, symab, symp, symq

rc = 0

!T test, if there are any ab pair

if (nvb(syma)*nvb(symb) == 0) then
  rc = 1
  ! RC=1 : there are no ab pair in this symmetry
  return
end if

! def initial address

irec = 1
symab = mul(syma,symb)

! loop over all combinations

do symp=1,nsym
  symq = mul(symab,symp)

  ! define number of records, required to store this block
  ! and determine shift in initial positions

  lengthpq = norb(symp)*norb(symq)
  nrecc = int(lengthpq/reclen)
  nrest = lengthpq-nrecc*reclen
  if (nrest > 0) nrecc = nrecc+1

  do a=1,nvb(syma)

    if (syma == symb) then
      bup = a
    else
      bup = nvb(symb)
    end if

    do b=1,bup

      abmap(a,b,symp) = irec
      irec = irec+nrecc

    end do
  end do
end do

return

end subroutine mkabpqmap
