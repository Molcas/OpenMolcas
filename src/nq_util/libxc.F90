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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************

module libxc

use Definitions, only: wp, iwp

implicit none
private

real(kind=wp), allocatable :: dfunc_dlapl(:,:), dfunc_drho(:,:), dfunc_dsigma(:,:), dfunc_dtau(:,:), func(:)
logical(kind=iwp) :: Only_exc = .false.

public :: dfunc_dlapl, dfunc_drho, dfunc_dsigma, dfunc_dtau, func, Only_exc

end module libxc
