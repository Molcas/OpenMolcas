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

subroutine set0(wrk,wrksize,map)
! this routine vanishes given mediate

use ccsd_global, only: Map_Type
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: map
integer(kind=iwp) :: ii, pos0

!1 def pos0
pos0 = map%d(1,1)

!2 set appropriate wrk = 0
ii = map%d(0,5)
wrk(pos0:map%d(ii,1)+map%d(ii,2)-1) = Zero

return

end subroutine set0
