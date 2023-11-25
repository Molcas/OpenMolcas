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

subroutine grad_cvb(grad)

use casvb_global, only: civb2, civb3, civb4, civb6, civb7, civb8, dvbdet, dxmove, grad1, grad2, gradx, icrit, memplenty, ovraa, &
                        ovraa_try, ovrab, ovrab_try, vec1, ww, ww_try
use Definitions, only: wp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: grad(*)

call touch_cvb('OOHESS')
if (dxmove .and. memplenty) then
  call cicopy_cvb(civb6,civb2)
  call cicopy_cvb(civb7,civb3)
  call cicopy_cvb(civb8,civb4)
else if (dxmove) then
  call cird_cvb(civb2,61006.2_wp)
  call cird_cvb(civb3,61007.2_wp)
  call cird_cvb(civb4,61008.2_wp)
end if
ovraa = ovraa_try
ovrab = ovrab_try
ww = ww_try
if (icrit == 1) then
  call gr_svb1_cvb(civb2,civb3,civb4,dvbdet,grad,grad1,grad2,gradx,vec1)
else if (icrit == 2) then
  call gr_evb1_cvb(civb2,civb3,civb4,dvbdet,grad,grad1,grad2,gradx,vec1)
end if

return

end subroutine grad_cvb
