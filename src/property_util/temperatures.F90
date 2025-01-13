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

module Temperatures

! Default temperatures for thermochemistry (MCLR, SLAPAF)

use Constants, only: Zero
use Definitions, only: wp

implicit none
private

real(kind=wp), parameter :: DefTemp(7) = [Zero,100.0_wp,273.15_wp,298.15_wp,323.15_wp,373.15_wp,473.15_wp]

public :: DefTemp

end module Temperatures
