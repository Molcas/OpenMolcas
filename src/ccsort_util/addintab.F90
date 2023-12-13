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

use ccsort_global, only: lunab, lunda1, map3, mbas, NORB, NSYM, nvb, reclen
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize, syma, symb, abmap(mbas,mbas,8)
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
integer(kind=iwp) :: a, b, bup, ii, irec0, length, pos, pos3, rc, symab, symp, symq

! def symab
symab = mul(syma,symb)

! make map3 for <_a_b|pq>

! set map3%i=0 (partly)

map3%i(1:nsym,1:nsym,1:nsym) = 0

! def 0-th row

map3%d(0,1) = 5
map3%d(0,2) = 5
map3%d(0,3) = 0
map3%d(0,4) = 0
map3%d(0,5) = nsym
map3%d(0,6) = 0

! def other rows

pos = map3%pos0
do ii=1,nsym

  symp = ii
  symq = mul(symab,symp)
  length = norb(symp)*norb(symq)
  map3%d(ii,1) = pos
  map3%d(ii,2) = length
  map3%d(ii,3) = symp
  map3%d(ii,4) = symq
  map3%d(ii,5) = 1
  map3%d(ii,6) = 1
  map3%i(symp,1,1) = ii
  pos = pos+length

end do

! write mapd,mapi to INTAB
call dawrtmap(lunab,map3,rc)

!T if there are no _a_b,pq integrals in this symab, skip summation over ab

if ((map3%d(nsym,1)+map3%d(nsym,2)) == map3%pos0) return

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
      ii = map3%i(symp,1,1)
      pos3 = map3%d(ii,1)
      length = map3%d(ii,2)

      ! read this block to #3
      if (length > 0) call daread(lunda1,irec0,wrk(pos3),length,reclen)

    end do

    ! since there must be some integrals, write them to TEMPAB

    call deflength(map3,length)
    call dawri(lunab,length,wrk(map3%pos0))

  end do
end do

return

end subroutine addintab
