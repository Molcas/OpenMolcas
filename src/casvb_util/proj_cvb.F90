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

!*********************************************************************
!*                                                                   *
!*  PROJ      := Project CASSCF vector onto irrep(s).                *
!*                                                                   *
!*********************************************************************
subroutine proj_cvb(civec)

use casvb_global, only: mxirrep, ndet, projsym
use Definitions, only: wp

implicit none
real(kind=wp), intent(inout) :: civec(0:ndet)
real(kind=wp) :: dum(mxirrep)

if (projsym) call psym1_cvb(civec(1:),civec(1:),dum,1)

return

end subroutine proj_cvb
