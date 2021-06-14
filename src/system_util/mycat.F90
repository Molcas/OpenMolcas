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

subroutine mycat(Junk,str0,str1,str2)
! Junk=str0//str1//str2

character*(*) Junk, str0, str1, str2

maxlen = len(Junk)

Junk = ' '
ile = 1
il = mylen(str0)
if (il > 0) then
  ils = 1
  ile = il+1
  if (ile > maxlen) goto 100
  Junk(ils:ile) = str0(1:il)
end if
il = mylen(str1)
if (il > 0) then
  ils = ile+1
  ile = ile+il
  if (ile > maxlen) goto 100
  Junk(ils:ile) = str1(1:il)
end if
il = mylen(str2)
if (il > 0) then
  ils = ile+1
  ile = ile+il
  if (ile > maxlen) goto 100
  Junk(ils:ile) = str2(1:il)
end if
return
100 write(6,*) ' too long strings to concatenate: '
write(6,*) str0,str1,str2

return

end subroutine mycat
