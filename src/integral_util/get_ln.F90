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

function get_ln(lunit)

use getline_mod, only: Quit_On_Error
use Definitions, only: iwp

implicit none
character(len=180) :: get_ln
integer(kind=iwp), intent(in) :: lunit
character(len=180), external :: get_ln_quit

get_ln = get_ln_quit(lunit,1)
if (Quit_On_Error) then
  call WarningMessage(2,'Error in Get_Ln')
  call Quit_OnUserError()
end if

return

end function get_ln
