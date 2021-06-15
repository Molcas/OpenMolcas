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

function mylen(s)
! return real length of the string without spaces...

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: mylen
character(len=*), intent(in) :: s
integer(kind=iwp) :: i, il

il = len(s)
if (il == 0) then
  mylen = 0
  return
end if
do i=il,1,-1
  if (s(i:i) /= ' ') then
    mylen = i
    return
  end if
end do
mylen = 0

return

end function mylen
