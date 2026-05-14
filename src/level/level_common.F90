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
! Copyright (C) 2022, Nike Dattani                                     *
!***********************************************************************

module level_common

use Definitions, only: wp, iwp

implicit none
private

real(kind=wp) :: ARV, PRV
real(kind=wp), allocatable :: DRDY2(:), RVB(:), SDRDY(:), VBZ(:), YVB(:)
integer(kind=iwp), parameter :: NDIMR = 200001
! Damping function parameters for printout .....
! what are these numbers?
real(kind=wp), parameter :: bDS(-4:0) = [2.50_wp,2.90_wp,3.3_wp,3.69_wp,3.95_wp], &
                            bTT(-1:2) = [2.44_wp,2.78_wp,3.126_wp,3.471_wp], &
                            cDS(-4:0) = [0.468_wp,0.446_wp,0.423_wp,0.40_wp,0.39_wp]

public :: ARV, bDS, bTT, cDS, DRDY2, NDIMR, PRV, RVB, SDRDY, VBZ, YVB

end module level_common
