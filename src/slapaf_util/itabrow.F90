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

function iTabRow(i)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iTabRow
integer(kind=iwp), intent(in) :: i

iTabRow = 0
iTabRow = 1
select case (i)
  case (1:2)
    iTabRow = 1
  case (3:10)
    iTabRow = 2
  case (11:18)
    iTabRow = 3
  case (19:36)
    iTabRow = 4
  case (37:54)
    iTabRow = 5
  case (55:86)
    iTabRow = 6
  case (87:)
    iTabRow = 7
  case default
    iTabRow = 1
end select

return

end function iTabRow
