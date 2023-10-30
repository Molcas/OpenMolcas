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

!***********************************************************************
!*                                                                     *
!*  SYMWEIGHT := CASSCF scalar product divided into irreps.            *
!*                                                                     *
!***********************************************************************
subroutine symweight_cvb(civec1,civec2,osym)

use casvb_global, only: mxirrep, ndet
use Definitions, only: wp

implicit none
real(kind=wp), intent(inout) :: civec1(0:ndet)
real(kind=wp), intent(in) :: civec2(0:ndet)
real(kind=wp), intent(out) :: osym(mxirrep)

call psym1_cvb(civec1(1:),civec2(1:),osym,2)

return

end subroutine symweight_cvb
