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

implicit none
real*8 DIASH(*)
real*8 SMAX
integer ISHLAB
#include "cholesky.fh"
integer JSHLAB
character*13 SECNAM
parameter(SECNAM='CHO_GETMAXSHL')

SMAX = -1.0d9
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
  DIASH(ISHLAB) = 0.0d0
end if

end subroutine CHO_GETMAXSHL
