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

function iPL_espf()
! Returns the print level

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iPL_espf
integer(kind=iwp) :: iPL
integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: Reduce_Prt

iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = 0
iPL_espf = iPL

return

end function iPL_espf
