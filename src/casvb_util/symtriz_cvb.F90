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

subroutine symtriz_cvb(orbs,cvb)

use casvb_global, only: norb, nvb
use Definitions, only: wp

implicit none
real(kind=wp), intent(inout) :: orbs(norb,norb), cvb(nvb)

call symtrizorbs_cvb(orbs)
call symtrizcvb_cvb(cvb)

return

end subroutine symtriz_cvb
