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

#include "reorg.fh"
#include "ccsort.fh"

integer abmap(1:mbas,1:mbas,1:8)
integer syma, symb, rc
! help variables
integer a, b, bup, symp, symq, symab
integer lengthpq, nrecc, nrest, irec

!T test, if there are any ab pair

if (nvb(syma)*nvb(symb) == 0) then
  rc = 1
  ! RC=1 : there are no ab pair in this symmetry
  return
else
  rc = 0
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
  nrecc = int(lengthpq/recl)
  nrest = lengthpq-nrecc*recl
  if (nrest > 0) then
    nrecc = nrecc+1
  end if

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
