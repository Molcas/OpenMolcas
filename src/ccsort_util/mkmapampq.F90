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

#include "ccsort.fh"
#include "reorg.fh"
integer syma
! help variables
integer symm, symp, symq, symmp
integer nhelp, position, length

! set mapi1 to zero

do symq=1,nsym
  do symp=1,nsym
    do symm=1,nsym
      mapi2(symm,symp,symq) = 0
    end do
  end do
end do

! def zero-th row

mapd2(0,1) = 1
mapd2(0,2) = 5
mapd2(0,3) = 5
mapd2(0,4) = 0
mapd2(0,6) = 0

nhelp = 0
position = poss20
do symm=1,nsym
  do symp=1,nsym
    symmp = mul(symm,symp)
    symq = mul(syma,symmp)
    nhelp = nhelp+1

    ! calc. length
    length = noa(symm)*NORB(symp)*NORB(symq)

    mapd2(nhelp,1) = position
    mapd2(nhelp,2) = length
    mapd2(nhelp,3) = symm
    mapd2(nhelp,4) = symp
    mapd2(nhelp,5) = symq
    mapd2(nhelp,6) = 1
    position = position+length

    mapi2(symm,symp,1) = nhelp

  end do
end do

mapd2(0,5) = nhelp

return

end subroutine mkmapampq
