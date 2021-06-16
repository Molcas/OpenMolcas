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
! Copyright (C) 2001, Valera Veryazov                                  *
!***********************************************************************
!***********************************************************************
!  SysWarnMsg
!
!> @brief
!>   Print nice formatted warning message
!> @author V. Veryazov
!>
!> @details
!> Print formatted message.
!> For a set of standard messages (started from ``MSG:``) the aliases
!> (defined in ::SysExpand) will be used.
!>
!> @param[in] Location routine name
!> @param[in] Text1    message text
!> @param[in] Text2    message text
!***********************************************************************

subroutine SysWarnMsg(Location,Text1,Text2)

use warnings, only: MaxWarnMess
use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: Location, Text1, Text2
integer(kind=iwp) :: i, Level
character(len=256) :: Str

! for these messages assume that Level is 1
Level = 1
if (Level > MaxWarnMess) MaxWarnMess = Level
call SysPutsStart()
call SysPuts('Location: ',Location,'\n\n\n')
call SysExpand(Text1,Str,i)
if (i == 0) then
  call SysPuts(Text1,' ',Text2)
else
  call SysPuts(str(:i),' ',Text2)
end if
call SysPutsEnd()

return

end subroutine SysWarnMsg
