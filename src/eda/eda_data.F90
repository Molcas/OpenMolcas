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
! Copyright (C) 2026, Okan Koeksal                                     *
!***********************************************************************

module EDA_Data

use Definitions, only: iwp

implicit none
private

integer(kind=iwp), public :: EDA_nFrag = 0
integer(kind=iwp), allocatable, public :: EDA_nAtom(:)
integer(kind=iwp), allocatable, public :: EDA_Charge(:)
integer(kind=iwp), allocatable, public :: EDA_Mult(:)
character(len=8), public :: EDA_RefMethod = 'AUTO    '
character(len=8), public :: EDA_FragBasis = 'OWN     '
logical, public :: EDA_Debug = .false.

public :: EDA_Reset_Input, EDA_Allocate_Fragments

contains

subroutine EDA_Reset_Input()

  if (allocated(EDA_nAtom)) deallocate(EDA_nAtom)
  if (allocated(EDA_Charge)) deallocate(EDA_Charge)
  if (allocated(EDA_Mult)) deallocate(EDA_Mult)

  EDA_nFrag = 0
  EDA_RefMethod = 'AUTO    '
  EDA_FragBasis = 'OWN     '
  EDA_Debug = .false.

end subroutine EDA_Reset_Input

subroutine EDA_Allocate_Fragments(nFrag,ierr)

  integer(kind=iwp), intent(in) :: nFrag
  integer(kind=iwp), intent(out) :: ierr
  integer :: AllocStat

  ierr = 1
  if (nFrag < 1) return

  if (allocated(EDA_nAtom)) deallocate(EDA_nAtom)
  if (allocated(EDA_Charge)) deallocate(EDA_Charge)
  if (allocated(EDA_Mult)) deallocate(EDA_Mult)

  allocate(EDA_nAtom(nFrag),EDA_Charge(nFrag),EDA_Mult(nFrag),stat=AllocStat)
  if (AllocStat /= 0) return

  EDA_nAtom(:) = 0
  EDA_Charge(:) = 0
  EDA_Mult(:) = 1
  ierr = 0

end subroutine EDA_Allocate_Fragments

end module EDA_Data
