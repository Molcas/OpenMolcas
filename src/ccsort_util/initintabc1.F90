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

subroutine initintabc1()
! this routine writes corresponding mapd and mapi to INTAB
! for nonsymmetrical (C1) case

use ccsort_global, only: lunab, map3, NORB, NSYM
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ii, length, pos, rc, syma, symab, symb, symp, symq

! def symab
syma = 1
symb = 1
symab = mul(syma,symb)

! make %d,%i for <_a_b|pq>

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

return

end subroutine initintabc1
