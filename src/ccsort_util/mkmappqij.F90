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
! this routine prepares mapd,mapi
! for <pq|ij> for p,q, i>=j to mapd1,mapi1

use ccsort_global, only: mapd1, mapi1, noa, NORB, NSYM, pos10
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: length, nhelp, pos, symi, symj, symp, sympq, sympqi, symq

! set mapi1 to zero

mapi1(1:nsym,1:nsym,1:nsym) = 0

! def zero-th row

mapd1(0,1) = 5
mapd1(0,2) = 5
mapd1(0,3) = 1
mapd1(0,4) = 1
mapd1(0,6) = 3

nhelp = 0
pos = pos10
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

      mapd1(nhelp,1) = pos
      mapd1(nhelp,2) = length
      mapd1(nhelp,3) = symp
      mapd1(nhelp,4) = symq
      mapd1(nhelp,5) = symi
      mapd1(nhelp,6) = symj
      pos = pos+length

      mapi1(symp,symq,symi) = nhelp

    end do
  end do
end do

mapd1(0,5) = nhelp

return

end subroutine mkmappqij
