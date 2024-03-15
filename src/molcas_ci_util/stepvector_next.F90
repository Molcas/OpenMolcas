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

subroutine STEPVECTOR_NEXT(MV,IDWN,IUP,STEPVECTOR)

use Definitions, only: iwp, u6
use gugx, only: NLEV, IOW1, NOW1

implicit none
integer(kind=iwp), intent(inout) :: MV, IDWN, IUP
integer(kind=iwp), intent(out) :: STEPVECTOR(NLEV)

! stop when MV is zero
if (MV == 0) then
  write(u6,'(1X,A)') 'stepvector_next has been depleted'
end if

call GETSTEPVECTOR(NOW1,IOW1,MV,IDWN,IUP,STEPVECTOR)

end subroutine STEPVECTOR_NEXT
