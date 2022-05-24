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

module Quater_globals

use Definitions, only: wp, iwp

implicit none
private

! debug
logical(kind=iwp) :: debug
! options
logical(kind=iwp) :: translate, rotate
! geoms
integer(kind=iwp), parameter ::MAXGEOMS=100
integer(kind=iwp) :: ngeoms, XYZ1(3), XYZ2(3)
type geoitem
  integer(kind=iwp) :: nat
  real(kind=wp), allocatable :: geo(:,:)
  character(len=20) :: title
  character(len=20), allocatable :: geolbl(:)
end type geoitem
type(geoitem) :: list(MAXGEOMS)
! rotation
logical(kind=iwp) :: matrixSet
real(kind=wp) :: RotMatrix(3,3)

public :: debug
public :: translate, rotate
public :: ngeoms, list, XYZ1, XYZ2
public :: matrixSet, RotMatrix

end module Quater_globals
