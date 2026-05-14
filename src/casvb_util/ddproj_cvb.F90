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

subroutine ddproj_cvb(c,nparm1)

use casvb_global, only: imethod, nprvb
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: c(*)
integer(kind=iwp), intent(in) :: nparm1
integer(kind=iwp) :: ioffs

if (.not. ((imethod == 4) .or. (imethod == 12))) call orthcvb_cvb(c,nparm1)
ioffs = nparm1-nprvb+1
call symtrizcvb_cvb(c(ioffs))

return

end subroutine ddproj_cvb
