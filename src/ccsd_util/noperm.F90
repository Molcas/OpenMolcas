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

subroutine noperm(wrk,wrksize,mapda,mapia,mapdb,mapib,poss0,posst)
! realize mapping without permutation
! define mapd,mapi

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, mapda(0:512,6), mapia(8,8,8), mapdb(0:512,6), mapib(8,8,8), poss0, posst
real(kind=wp) :: wrk(wrksize)
#include "ccsd1.fh"
integer(kind=iwp) :: i, ib, j, k, nhelp

! def mapib

do k=1,nsym
  do j=1,nsym
    do i=1,nsym
      mapib(i,j,k) = mapia(i,j,k)
    end do
  end do
end do

! def initial values

do nhelp=1,6
  mapdb(0,nhelp) = mapda(0,nhelp)
end do

posst = poss0
do ib=1,mapda(0,5)
  do nhelp=2,6
    mapdb(ib,nhelp) = mapda(ib,nhelp)
  end do
  mapdb(ib,1) = posst
  posst = posst+mapdb(ib,2)

  call map11(wrk(mapda(ib,1)),wrk(mapdb(ib,1)),mapda(ib,2),1)

end do

return

end subroutine noperm
