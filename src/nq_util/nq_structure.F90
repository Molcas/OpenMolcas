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

module NQ_Structure

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

type Info_Ang_type
  integer(kind=iwp) :: L_eff = 0
  integer(kind=iwp) :: nPoints = 0
  real(kind=wp), allocatable :: R(:,:)
end type Info_Ang_type

type NQ_data_type
  real(kind=wp), allocatable :: Coor(:)
  real(kind=wp) :: A_High = -huge(Zero)
  real(kind=wp) :: A_Low = huge(Zero)
  real(kind=wp) :: R_RS = Zero
  real(kind=wp) :: R_max = Zero
  integer(kind=iwp) :: l_max = -1
  real(kind=wp), allocatable :: R_Quad(:,:)
  integer(kind=iwp), allocatable :: Angular(:)
  integer(kind=iwp) :: Atom_Nr = -1
  real(kind=wp), allocatable :: dOdx(:,:,:)
end type NQ_data_type

integer(kind=iwp), parameter :: LMax_NQ = 62
type(Info_Ang_type) Info_Ang(LMax_NQ)
type(NQ_data_type), allocatable :: NQ_data(:)

public :: Close_Info_Ang, Close_NQ_Data, Info_Ang, LMax_NQ, NQ_data

contains

subroutine Close_Info_Ang()

  use stdalloc, only: mma_deallocate

  integer(kind=iwp) iAngular

  do iAngular=1,size(Info_Ang)
    Info_Ang(iAngular)%L_eff = 0
    Info_Ang(iAngular)%nPoints = 0
    if (allocated(Info_Ang(iAngular)%R)) call mma_deallocate(Info_Ang(iAngular)%R)
  end do

end subroutine Close_Info_Ang

subroutine Close_NQ_Data()

  use stdalloc, only: mma_deallocate

  integer(kind=iwp) iNQ

  ! Cleanup and close
  do iNQ=1,size(NQ_data)
    call mma_deallocate(NQ_data(iNQ)%Coor)
    NQ_data(iNQ)%A_High = -huge(Zero)
    NQ_data(iNQ)%A_Low = huge(Zero)
    NQ_data(iNQ)%R_RS = Zero
    NQ_data(iNQ)%R_max = Zero
    NQ_data(iNQ)%l_Max = -1
    if (allocated(NQ_data(iNQ)%R_Quad)) call mma_deallocate(NQ_data(iNQ)%R_Quad)
    if (allocated(NQ_data(iNQ)%Angular)) call mma_deallocate(NQ_data(iNQ)%Angular)
    NQ_Data(iNQ)%Atom_Nr = -1
    if (allocated(NQ_data(iNQ)%dOdx)) call mma_deallocate(NQ_data(iNQ)%dOdx)
  end do
  deallocate(NQ_Data) !IFG

end subroutine Close_NQ_Data

end module NQ_Structure
