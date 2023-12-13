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
! Copyright (C) 2010, Francesco Aquilante                              *
!***********************************************************************

function Vt_lim(rho,drho,ddrho)

use Constants, only: One, Two, Quart
use Definitions, only: wp

implicit none
real(kind=wp) :: Vt_lim
real(kind=wp), intent(in) :: rho, drho(3), ddrho
real(kind=wp) :: rhoinv, rhoinv2, xnorm

rhoinv = One/rho
rhoinv2 = rhoinv**Two
xnorm = drho(1)**2+drho(2)**2+drho(3)**2

Vt_lim = 0.125_wp*xnorm*rhoinv2-Quart*ddrho*rhoinv

end function Vt_lim
