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

subroutine stz(wrk,wrksize,mapda)
! this routine vanishes A
! A = 0
!
! mapda - direct map of A m(I/O)
!
! N.B. this routine should be done using matrix operations

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, mapda(0:512,6)
real(kind=wp) :: wrk(wrksize)
integer(kind=iwp) :: nhelp1, nhelp2, nhelp3

!1 def the length of the mediate
nhelp1 = mapda(0,5)
nhelp3 = mapda(nhelp1,1)+mapda(nhelp1,2)-mapda(1,1)

!2 def initial position
nhelp2 = mapda(1,1)

!3 refactoring
do nhelp1=nhelp2,nhelp2+nhelp3-1
  wrk(nhelp1) = Zero
end do

return

end subroutine stz
