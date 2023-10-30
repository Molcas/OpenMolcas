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

subroutine applythmes_cvb(civbh,orbs)
! Apply T(O) (H - E) T(O) to CIVBH:

use casvb_global, only: gjorb, gjorb2, ndet, norb, ovraa, ww
use Definitions, only: wp

implicit none
real(kind=wp), intent(inout) :: civbh(0:ndet)
real(kind=wp), intent(in) :: orbs(norb,norb)

call makegjorbs_cvb(orbs)

call applyt_cvb(civbh,gjorb)
call proj_cvb(civbh)
call applyhpcx_cvb(civbh,-ww/ovraa)
call applyt_cvb(civbh,gjorb2)

return

end subroutine applythmes_cvb
