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

subroutine cct3_noperm(wrk,wrksize,mapda,mapia,mapdb,mapib,pos0,post)
! realize mapping without permutation
! define mapd,mapi

use CCT3_global, only: nsym
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, mapda(0:512,6), mapia(8,8,8), mapdb(0:512,6), mapib(8,8,8), pos0, post
real(kind=wp) :: wrk(wrksize)
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

post = pos0
do ib=1,mapda(0,5)
  do nhelp=2,6
    mapdb(ib,nhelp) = mapda(ib,nhelp)
  end do
  mapdb(ib,1) = post
  post = post+mapdb(ib,2)

  call cct3_map11(wrk(mapda(ib,1)),wrk(mapdb(ib,1)),mapda(ib,2),1)

end do

return

end subroutine cct3_noperm
