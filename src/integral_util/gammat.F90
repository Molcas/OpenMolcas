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

function gammat(x)
!***********************************************************************
!                                                                      *
! Object: to compute the angular contribution to the multipole integral*
!         between continuum basis functions within an R-matrix sphere  *
!         (theta integration)                                          *
!                                                                      *
!***********************************************************************

use rmat, only: m_Gam, n_Gam
use Constants, only: Zero, One, Two, Three, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: gammat
real(kind=wp), intent(in) :: x
integer(kind=iwp) :: lSinT, lCosT, k
real(kind=wp) :: arg1, arg2, arg3
real(kind=wp), external :: dGamma_Molcas

lsint = m_gam
lcost = n_gam
k = (-1)**lcost
if (k == (-1)) then
  gammat = Zero
else
  arg1 = Half*(real(lcost,kind=wp)+One)
  arg2 = Half*(real(lsint,kind=wp)+Two)
  arg3 = Half*(real(lsint,kind=wp)+real(lcost,kind=wp)+Three)
  gammat = dgamma_molcas(arg1)*dgamma_molcas(arg2)/dgamma_molcas(arg3)
end if

return
! Avoid unused argument warnings
if (.false.) call Unused_real(x)

end function gammat
