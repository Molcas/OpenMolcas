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

subroutine CHO_GETMAXSHL(DIASH,SMAX,ISHLAB)
!
! Purpose: Get max. shell pair and update DIASH.

use Cholesky, only: nnShl
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: DIASH(*)
real(kind=wp), intent(out) :: SMAX
integer(kind=iwp), intent(out) :: ISHLAB
integer(kind=iwp) :: JSHLAB
character(len=*), parameter :: SECNAM = 'CHO_GETMAXSHL'

SMAX = -1.0e9_wp
ISHLAB = -1
do JSHLAB=1,NNSHL
  if (DIASH(JSHLAB) > SMAX) then
    SMAX = DIASH(JSHLAB)
    ISHLAB = JSHLAB
  end if
end do

if (ISHLAB < 1) then
  call CHO_QUIT('Error in '//SECNAM,104)
else
  DIASH(ISHLAB) = Zero
end if

end subroutine CHO_GETMAXSHL
