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
! this routine prepares mapd,mapi
! for <am|pq> for given syma, m, p,q to mapd2,mapi2

use ccsort_global, only: mapd2, mapi2, noa, NORB, NSYM, pos20
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: syma
integer(kind=iwp) :: length, nhelp, pos, symm, symmp, symp, symq

! set mapi1 to zero

mapi2(1:nsym,1:nsym,1:nsym) = 0

! def zero-th row

mapd2(0,1) = 1
mapd2(0,2) = 5
mapd2(0,3) = 5
mapd2(0,4) = 0
mapd2(0,6) = 0

nhelp = 0
pos = pos20
do symm=1,nsym
  do symp=1,nsym
    symmp = mul(symm,symp)
    symq = mul(syma,symmp)
    nhelp = nhelp+1

    ! calc. length
    length = noa(symm)*NORB(symp)*NORB(symq)

    mapd2(nhelp,1) = pos
    mapd2(nhelp,2) = length
    mapd2(nhelp,3) = symm
    mapd2(nhelp,4) = symp
    mapd2(nhelp,5) = symq
    mapd2(nhelp,6) = 1
    pos = pos+length

    mapi2(symm,symp,1) = nhelp

  end do
end do

mapd2(0,5) = nhelp

return

end subroutine mkmapampq
