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

subroutine SysFileMsg(Location,Text1,Lu,Text2)

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: Location, Text1, Text2
integer(kind=iwp), intent(in) :: Lu
integer(kind=iwp) :: i
character(len=256) :: str

call SysPutsStart()
call SysPuts('Location: ',Location,'\n')
write(str,*) Lu
call SysPuts('Unit    : ',str,' ')
str = ' '
inquire(unit=lu,name=str)
if (str /= ' ') then
  call SysPuts('File    : ',str,'\n')
end if

call SysExpand(Text1,Str,i)
if (i == 0) then
  call SysPuts(Text1,'\n',Text2)
else
  call SysPuts(str(:i),'\n',Text2)
end if
call SysPutsEnd()
call Abend()

return

end subroutine SysFileMsg
