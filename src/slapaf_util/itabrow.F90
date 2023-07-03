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

integer function iTabRow(i)

implicit integer(a-z)

iTabRow = 0
iTabRow = 1
if ((i > 0) .and. (i <= 2)) then
  iTabRow = 1
else if ((i > 2) .and. (i <= 10)) then
  iTabRow = 2
else if ((i > 10) .and. (i <= 18)) then
  iTabRow = 3
else if ((i > 18) .and. (i <= 36)) then
  iTabRow = 4
else if ((i > 36) .and. (i <= 54)) then
  iTabRow = 5
else if ((i > 54) .and. (i <= 86)) then
  iTabRow = 6
else if (i > 86) then
  iTabRow = 7
end if

return

end function iTabRow
