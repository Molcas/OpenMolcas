!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

function find_Lu(filename)
! find the unit number assigned to an opened file (not translated)
! returns -1 if not found

use Fast_IO, only: isOpen, LuName, MxFile
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: find_Lu
character(len=*), intent(in) :: filename
integer(kind=iwp) :: n

find_Lu = -1
do n=1,MxFile
  if (isOpen(n) == 0) cycle
  if (LuName(n) == filename) then
    find_Lu = n
    exit
  end if
end do

return

end function find_Lu
