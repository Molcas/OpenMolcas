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

subroutine calcrh1(wrk,wrksize,v1,v2)
! this routine calcs V1 = V1-V2
!
! v1 - map type of vector 1 (I)
! v2 - map type of vector 2 (I)
!
! N.B. it is assumed, that V1 and V2 are of the same type

use ccsd_global, only: Map_Type
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: v1, v2
integer(kind=iwp) :: ii, length, pos1, pos2

!1 calc length
ii = v1%d(0,5)
length = v1%d(ii,1)+v1%d(ii,2)-v1%d(1,1)

!2 realize substract
if (length > 0) then
  pos1 = v1%d(1,1)
  pos2 = v2%d(1,1)
  wrk(pos1:pos1+length-1) = wrk(pos1:pos1+length-1)-wrk(pos2:pos2+length-1)
end if

return

end subroutine calcrh1
