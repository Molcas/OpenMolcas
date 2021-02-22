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
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************

subroutine Dispersion_Kriging(x0_,y_,ndimx)

use kriging_mod, only: sigma, x0
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ndimx
real(kind=wp), intent(in) :: x0_(ndimx)
real(kind=wp), intent(out) :: y_
real(kind=wp), parameter :: invnormal = 1.9599639845400542355_wp

!x0 is the n-dimensional vector of the coordinates for which the dispersion is computed

x0(:) = x0_(:)
call covarvector(0) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
call predict(0)
! 95% confidence -> 1.96*sigma
y_ = invnormal*sigma

return

end subroutine Dispersion_Kriging
