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

use ccsort_global, only: lunab, mapd3, mapi3, NORB, NSYM, pos30
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ii, length, pos, rc, syma, symab, symb, symp, symq

! def symab
syma = 1
symb = 1
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

return

end subroutine initintabc1
