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

subroutine SysWarnFileMsg(Location,TheFile,Text1,Text2)

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: Location, Text1, Text2, TheFile
integer(kind=iwp) :: i
character(len=256) :: Str

call SysPutsStart()
call SysPuts('Location: ',Location,'\n')
call SysExpand(TheFile,Str,i)
call SysPuts('File: ',TheFile,'\n\n\n')
call SysExpand(Text1,Str,i)
if (i == 0) then
  call SysPuts(Text1,' ',Text2)
else
  call SysPuts(str(:i),' ',Text2)
end if
call SysPutsEnd()

return

end subroutine SysWarnFileMsg
