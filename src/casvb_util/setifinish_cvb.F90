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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine setifinish_cvb(icode)

use casvb_global, only: ifinish
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: icode

if (icode == 1) then
  ifinish = 0
else if (icode == 3) then
  ifinish = 1
end if

return

end subroutine setifinish_cvb
