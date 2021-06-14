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

character*(*) str
character fmt*20

iTooLong = 60
i = len(str)
if (i > iTooLong+8) then
  ! oops! too long
  write(6,'(a,a)') ' ###    ',str
  return
end if
i = iTooLong+8-i
if (i == 0) then
  fmt = '(a,a,a)'
else
  write(fmt,'(a, i2,a)') '(a,a,',i,'x,a)'
end if
write(6,fmt) ' ###    ',str,' ###'

return

end subroutine SysDumpStr
