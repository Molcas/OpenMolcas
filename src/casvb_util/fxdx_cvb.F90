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

subroutine fxdx_cvb(fx,fast,dx)

use casvb_global, only: civb1, civb2, civb3, civb4, civb6, civb7, civb8, cvbdet, cvbtry, dxmove, icnt_ci, icrit, memplenty, orbstry
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: fx
logical(kind=iwp), intent(in) :: fast
real(kind=wp), intent(in) :: dx(*)

dxmove = .false.
call upd_cvb(dx,orbstry,cvbtry)
if (.not. memplenty) then
  call ciwr_cvb(civb2,61002.2_wp)
  call ciwr_cvb(civb3,61003.2_wp)
  call ciwr_cvb(civb4,61004.2_wp)
  icnt_ci(2:4) = 0
end if
icnt_ci(6:8) = 0
if (icrit == 1) then
  call fx_svb1_cvb(fx,fast,orbstry,cvbtry,civb1,civb6,civb7,civb8,cvbdet)
else if (icrit == 2) then
  call fx_evb1_cvb(fx,fast,orbstry,cvbtry,civb1,civb6,civb7,civb8,cvbdet)
end if
if (.not. memplenty) then
  call ciwr_cvb(civb6,61006.2_wp)
  call ciwr_cvb(civb7,61007.2_wp)
  call ciwr_cvb(civb8,61008.2_wp)
  call cird_cvb(civb2,61002.2_wp)
  call cird_cvb(civb3,61003.2_wp)
  call cird_cvb(civb4,61004.2_wp)
end if
! Figure out what we just calculated, and make it up2date:
if (dxmove) then
  if (icrit == 1) then
    call make_cvb('SVB')
  else if (icrit == 2) then
    call make_cvb('EVB')
  end if
else
  if (icrit == 1) then
    call make_cvb('SVBTRY')
  else if (icrit == 2) then
    call make_cvb('EVBTRY')
  end if
end if

return

end subroutine fxdx_cvb
