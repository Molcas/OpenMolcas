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

subroutine addintab(wrk,wrksize,syma,symb,abmap)
! this routine adds contributions to open INTAB1 file,
! comming from ab syma,symb

use ccsort_global, only: lunab, lunda1, mapd3, mapi3, mbas, NORB, NSYM, nvb, pos30, reclen
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize, syma, symb, abmap(mbas,mbas,8)
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
integer(kind=iwp) :: a, b, bup, ii, irec0, length, pos, pos3, rc, symab, symp, symq

! def symab
symab = mul(syma,symb)

! make mapd3,mapi3 for <_a_b|pq>

! set mapi3=0 (partly)

mapi3(1:nsym,1:nsym,1:nsym) = 0

! def 0-th row

mapd3(0,1) = 5
mapd3(0,2) = 5
mapd3(0,3) = 0
mapd3(0,4) = 0
mapd3(0,5) = nsym
mapd3(0,6) = 0

! def other rows

pos = pos30
do ii=1,nsym

  symp = ii
  symq = mul(symab,symp)
  length = norb(symp)*norb(symq)
  mapd3(ii,1) = pos
  mapd3(ii,2) = length
  mapd3(ii,3) = symp
  mapd3(ii,4) = symq
  mapd3(ii,5) = 1
  mapd3(ii,6) = 1
  mapi3(symp,1,1) = ii
  pos = pos+length

end do

! write mapd,mapi to INTAB
call dawrtmap(lunab,mapd3,mapi3,rc)

!T if there are no _a_b,pq integrals in this symab, skip summation over ab

if ((mapd3(nsym,1)+mapd3(nsym,2)) == pos30) return

! loop over a,b

do a=1,nvb(syma)

  if (syma == symb) then
    bup = a
  else
    bup = nvb(symb)
  end if

  do b=1,bup

    ! loop over symp

    do symp=1,nsym

      ! def irec0 for this a,b,symp in TEMPDA1
      irec0 = abmap(a,b,symp)

      ! def corresponding position and length in #3
      ii = mapi3(symp,1,1)
      pos3 = mapd3(ii,1)
      length = mapd3(ii,2)

      ! read this block to #3
      if (length > 0) call daread(lunda1,irec0,wrk(pos3),length,reclen)

    end do

    ! since there must be some integrals, write them to TEMPAB

    call deflength(mapd3,length)
    call dawri(lunab,length,wrk(pos30))

  end do
end do

return

end subroutine addintab
