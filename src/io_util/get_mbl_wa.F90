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

function get_MBl_wa()
! help function to get the value without using the module

use Fast_IO, only: MBl_wa
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: get_MBl_wa

get_MBl_wa = MBl_wa

return

end function get_MBl_wa
