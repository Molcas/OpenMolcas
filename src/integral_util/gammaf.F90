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

function gammaf(x)
!***********************************************************************
!                                                                      *
! Object: to compute the angular contribution to the multipole integral*
!         between continuum basis functions within an R-matrix sphere  *
!         (phi integration)                                            *
!                                                                      *
!***********************************************************************

use rmat, only: lcosf, lsinf, m_gam, n_gam
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: gammaf
real(kind=wp), intent(in) :: x
integer(kind=iwp) :: k1, k2
real(kind=wp) :: arg1, arg2, arg3
real(kind=wp), external :: dgamma_molcas

lcosf = n_gam
lsinf = m_gam
k1 = (-1)**lsinf
k2 = (-1)**lcosf
if ((k1 == -1) .or. (k2 == -1)) then
  gammaf = Zero
else
  arg1 = Half*(real(lcosf,kind=wp)+One)
  arg2 = Half*(real(lsinf,kind=wp)+One)
  arg3 = Half*(real(lsinf,kind=wp)+real(lcosf,kind=wp)+Two)
  gammaf = Two*dgamma_molcas(arg1)*dgamma_molcas(arg2)/dgamma_molcas(arg3)
end if

return
! Avoid unused argument warnings
if (.false.) call Unused_real(x)

end function gammaf
