!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine RPA_Warn(Level,Message)

! Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
! Level <= 1 --- issue warning with message Message and return.
! Level >= 2 --- issue warning with message Message and quit.
!                Error codes used in xQuit:
!                Level=2: _RC_INPUT_ERROR_
!                Level=3: _RC_INTERNAL_ERROR_
!                Level>3: _RC_GENERAL_ERROR_

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Level
character(len=*), intent(in) :: Message
#include "warnings.h"
integer(kind=iwp) :: iLevel, rc

if (Level <= 1) then
  iLevel = max(Level,0)
  rc = 0
else
  iLevel = 2
  if (Level == 2) then
    rc = _RC_INPUT_ERROR_
  else if (Level == 3) then
    rc = _RC_INTERNAL_ERROR_
  else
    rc = _RC_GENERAL_ERROR_
  end if
end if
call WarningMessage(iLevel,Message)
call xFlush(u6)
if (rc /= 0) then
  call xQuit(rc)
end if

end subroutine RPA_Warn
