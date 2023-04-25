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

use Definitions, only: wp

implicit none
private

real(kind=wp) :: ARV, PRV
real(kind=wp), allocatable :: DRDY2(:), RVB(:), SDRDY(:), VBZ(:), YVB(:)

public :: ARV, DRDY2, PRV, RVB, SDRDY, VBZ, YVB

end module level_common
