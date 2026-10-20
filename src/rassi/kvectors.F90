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

module kVectors

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: nk_Vector
real(kind=wp) :: e_Vector(3)
real(kind=wp), dimension(:,:), allocatable :: k_Vector

public :: e_Vector, k_Vector, nk_Vector

end module kVectors
