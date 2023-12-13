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

subroutine mkmapampq(syma)
! this routine prepares %d, %i
! for <am|pq> for given syma, m, p,q to map2

use ccsort_global, only: map2, noa, NORB, NSYM
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: syma
integer(kind=iwp) :: length, nhelp, pos, symm, symmp, symp, symq

! set %i to zero

map2%i(1:nsym,1:nsym,1:nsym) = 0

! def zero-th row

map2%d(0,1) = 1
map2%d(0,2) = 5
map2%d(0,3) = 5
map2%d(0,4) = 0
map2%d(0,6) = 0

nhelp = 0
pos = map2%pos0
do symm=1,nsym
  do symp=1,nsym
    symmp = mul(symm,symp)
    symq = mul(syma,symmp)
    nhelp = nhelp+1

    ! calc. length
    length = noa(symm)*NORB(symp)*NORB(symq)

    map2%d(nhelp,1) = pos
    map2%d(nhelp,2) = length
    map2%d(nhelp,3) = symm
    map2%d(nhelp,4) = symp
    map2%d(nhelp,5) = symq
    map2%d(nhelp,6) = 1
    pos = pos+length

    map2%i(symm,symp,1) = nhelp

  end do
end do

map2%d(0,5) = nhelp

return

end subroutine mkmapampq
