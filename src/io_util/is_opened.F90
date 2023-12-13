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

function is_opened(lu)
! help function to get the value without using the module

use Fast_IO, only: isOpen, MxFile
use Definitions, only: iwp

implicit none
logical(kind=iwp) :: is_opened
integer(kind=iwp), intent(in) :: lu

if ((lu > 0) .and. (lu < MxFile)) then
  is_opened = isOpen(lu) /= 0
else
  is_opened = .false.
end if

return

end function is_opened
