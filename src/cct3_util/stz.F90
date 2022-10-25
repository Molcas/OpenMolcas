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

subroutine stz(wrk,wrksize,a)
! this routine vanishes A
! A = 0
!
! a - A m(I/O)
!
! N.B. this routine should be done using matrix operations

use CCT3_global, only: Map_Type
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
type(Map_Type), intent(in) :: a
integer(kind=iwp) :: nhelp1, nhelp2, nhelp3

!1 def the length of the mediate
nhelp1 = a%d(0,5)
nhelp3 = a%d(nhelp1,1)+a%d(nhelp1,2)-a%d(1,1)

!2 def initial position
nhelp2 = a%d(1,1)

!3 refactoring
wrk(nhelp2:nhelp2+nhelp3-1) = Zero

return

end subroutine stz
