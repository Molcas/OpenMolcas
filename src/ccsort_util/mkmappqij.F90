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

#include "ccsort.fh"
#include "reorg.fh"
! help variables
integer symi, symj, symp, symq, sympq, sympqi
integer nhelp, position, length

! set mapi1 to zero

do symi=1,nsym
  do symq=1,nsym
    do symp=1,nsym
      mapi1(symp,symq,symi) = 0
    end do
  end do
end do

! def zero-th row

mapd1(0,1) = 5
mapd1(0,2) = 5
mapd1(0,3) = 1
mapd1(0,4) = 1
mapd1(0,6) = 3

nhelp = 0
position = poss10
do symp=1,nsym
  do symq=1,nsym
    sympq = mul(symp,symq)
    do symi=1,nsym
      sympqi = mul(sympq,symi)
      symj = sympqi
      if (symj > symi) goto 102
      nhelp = nhelp+1

      ! calc. length
      length = noa(symi)*noa(symj)*NORB(symp)*NORB(symq)

      mapd1(nhelp,1) = position
      mapd1(nhelp,2) = length
      mapd1(nhelp,3) = symp
      mapd1(nhelp,4) = symq
      mapd1(nhelp,5) = symi
      mapd1(nhelp,6) = symj
      position = position+length

      mapi1(symp,symq,symi) = nhelp

      102 continue
    end do
  end do
end do

mapd1(0,5) = nhelp

return

end subroutine mkmappqij
