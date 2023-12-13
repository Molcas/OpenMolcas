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

subroutine deflength(map,length)
! this routine defines length of mediate, described by map

use ccsort_global, only: Map_Type
use Definitions, only: iwp

implicit none
type(Map_Type), intent(in) :: map
integer(kind=iwp), intent(out) :: length
integer(kind=iwp) :: ii

ii = map%d(0,5)
length = map%d(ii,1)+map%d(ii,2)-map%d(1,1)

return

end subroutine deflength
