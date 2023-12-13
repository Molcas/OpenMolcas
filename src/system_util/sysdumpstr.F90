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

subroutine SysDumpStr(str)

use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: str
integer(kind=iwp) :: i, iTooLong
character(len=20) :: frmt

iTooLong = 60
i = len(str)
if (i > iTooLong+8) then
  ! oops! too long
  write(u6,'(a,a)') ' ###    ',str
  return
end if
i = iTooLong+8-i
if (i == 0) then
  frmt = '(a,a,a)'
else
  write(frmt,'(a, i2,a)') '(a,a,',i,'x,a)'
end if
write(u6,frmt) ' ###    ',str,' ###'

return

end subroutine SysDumpStr
