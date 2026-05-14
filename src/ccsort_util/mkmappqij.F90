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

subroutine mkmappqij()
! this routine prepares %d, %i
! for <pq|ij> for p,q, i>=j to map1

use ccsort_global, only: map1, noa, NORB, NSYM
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: length, nhelp, pos, symi, symj, symp, sympq, sympqi, symq

! set map1%i to zero

map1%i(1:nsym,1:nsym,1:nsym) = 0

! def zero-th row

map1%d(0,1) = 5
map1%d(0,2) = 5
map1%d(0,3) = 1
map1%d(0,4) = 1
map1%d(0,6) = 3

nhelp = 0
pos = map1%pos0
do symp=1,nsym
  do symq=1,nsym
    sympq = mul(symp,symq)
    do symi=1,nsym
      sympqi = mul(sympq,symi)
      symj = sympqi
      if (symj > symi) cycle
      nhelp = nhelp+1

      ! calc. length
      length = noa(symi)*noa(symj)*NORB(symp)*NORB(symq)

      map1%d(nhelp,1) = pos
      map1%d(nhelp,2) = length
      map1%d(nhelp,3) = symp
      map1%d(nhelp,4) = symq
      map1%d(nhelp,5) = symi
      map1%d(nhelp,6) = symj
      pos = pos+length

      map1%i(symp,symq,symi) = nhelp

    end do
  end do
end do

map1%d(0,5) = nhelp

return

end subroutine mkmappqij
