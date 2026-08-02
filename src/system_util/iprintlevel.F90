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
! Copyright (C) Valera Veryazov                                        *
!***********************************************************************
!  iPrintLevel
!
!> @brief
!>   Check or set global print level
!> @author V. Veryazov
!>
!> @details
!> If Level is a number &ge; ``0``, set print level,
!> else return current print level.
!>
!> In a first call of the routine environment variable
!> ``MOLCAS_PRINT`` is used to set the initial print level.
!>
!> Allowed values are:
!>
!> - ``SILENT`` (``0``)
!> - ``TERSE`` (``1``)
!> - ``NORMAL`` (``2``)
!> - ``VERBOSE`` (``3``)
!> - ``DEBUG`` (``4``)
!> - ``INSANE`` (``5``)
!>
!> @param[in] Level Print Level
!>
!> @return Print level
!***********************************************************************

function iPrintLevel(Level)

use PrintLevel, only: DEBUG, INSANE, SILENT, TERSE, USUAL, VERBOSE
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iPrintLevel
integer(kind=iwp), intent(in) :: Level
integer(kind=iwp), save :: inum, isFirst = 0, istatus, nPrintLevel
character(len=80) :: Val

if (Level >= 0) then
  nPrintLevel = Level
  isFirst = 1
  iPrintLevel = Level
  return
end if
if (isFirst == 0) then
  call getenvf('MOLCAS_PRINT',Val)
  call UpCase(Val)
  select case (Val)
    case ('SILENT')
      nPrintLevel = SILENT
    case ('TERSE')
      nPrintLevel = TERSE
    case ('NORMAL','USUAL')
      nPrintLevel = USUAL
    case ('VERBOSE')
      nPrintLevel = VERBOSE
    case ('DEBUG')
      nPrintLevel = DEBUG
    case ('INSANE')
      nPrintLevel = INSANE
    case default
      inum = -1
      read(Val,*,iostat=istatus) inum
      if ((istatus == 0) .and. (inum >= SILENT) .and. (inum <= INSANE)) then
        nPrintLevel = inum
      else
        nPrintLevel = USUAL
      end if
  end select
end if
iPrintLevel = nPrintLevel

return

end function iPrintLevel
