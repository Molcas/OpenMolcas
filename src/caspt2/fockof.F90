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

module Fockof

use definitions, only: iwp, wp

implicit none
private

public :: IOFFIT, IOFFIA, IOFFTA
public :: FIT, FIT_Full
public :: FTI, FTI_Full
public :: FIA, FIA_Full
public :: FAI, FAI_Full
public :: FTA, FTA_Full
public :: FAT, FAT_Full

type rPointers
  real(kind=wp), pointer :: A(:) => null()
end type rPointers

type(rPointers) :: FIT(8)
real(kind=wp), allocatable, target :: FIT_Full(:)
type(rPointers) :: FTI(8)
real(kind=wp), allocatable, target :: FTI_Full(:)
type(rPointers) :: FIA(8)
real(kind=wp), allocatable, target :: FIA_Full(:)
type(rPointers) :: FAI(8)
real(kind=wp), allocatable, target :: FAI_Full(:)
type(rPointers) :: FTA(8)
real(kind=wp), allocatable, target :: FTA_Full(:)
type(rPointers) :: FAT(8)
real(kind=wp), allocatable, target :: FAT_Full(:)

integer(kind=iwp) :: IOFFIT(8), IOFFIA(8), IOFFTA(8)

end module Fockof
