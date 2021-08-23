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
!***********************************************************************
!  WarningMessage
!
!> @brief
!>   Print warning message
!> @author V. Veryazov
!>
!> @details
!> Print message in uniform format.
!>
!> @note
!> Routine updates \c MaxWarnMess.
!>
!> Recommended levels:
!> - ``0``: for message without recording that it was an error
!> - ``1``: for warnings
!> - ``2``: for errors
!>
!> @param[in] Level Warning level
!> @param[in] STR   Message
!***********************************************************************

subroutine WarningMessage(Level,STR)

use warnings, only: MaxWarnMess
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: Level
character(len=*), intent(in) :: STR

if (Level > MaxWarnMess) MaxWarnMess = Level
call SysPutsStart()
if (Level == 1) then
  call SysPuts('WARNING: ',STR,' ')
else if (Level == 2) then
  call SysPuts('ERROR: ',STR,' ')
else
  call SysPuts(STR,' ',' ')
end if
call SysPutsEnd()

!write(u6,'(A)') '*** '
!jj = 1
!do
!  i = index(STR(jj:),';')
!  if (i == 0) exit
!  write(u6,'(A,A)') '*** ',STR(jj:jj+i-2)
!  jj = i+jj
!end do
!write(u6,'(A,A)') '*** ',STR(jj:)
!write(u6,'(A)') '*** '

return

end subroutine WarningMessage
