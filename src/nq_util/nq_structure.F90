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

implicit none
private
public :: NQ_data, Close_NQ_Data, Info_Ang, Close_Info_Ang, LMax_NQ

#include "stdalloc.fh"

!define declare_ip_dodx     ip_dOdx(iNQ,i) = ipNQ+(iNQ-1)*l_NQ+15+(iTabMx+1)+(i-1)*9

type NQ_data_raw
  sequence
  real*8, allocatable :: Coor(:)
  real*8 :: A_High = -1.0d99
  real*8 :: A_Low = 1.0d99
  real*8 :: R_RS = 0.0d0
  real*8 :: R_max = 0.0d0
  integer :: l_max = -1
  real*8, allocatable :: R_Quad(:,:)
  integer, allocatable :: Angular(:)
  integer :: Atom_Nr = -1
  real*8, allocatable :: dOdx(:,:,:)
end type NQ_data_raw

type(NQ_data_raw), allocatable :: NQ_data(:)

type Info_A
  sequence
  integer :: L_eff = 0
  integer :: nPoints = 0
  real*8, allocatable :: R(:,:)
end type Info_A

integer, parameter :: LMax_NQ = 62
type(Info_A) Info_Ang(LMax_NQ)

contains

subroutine Close_Info_Ang()

  integer iAngular

  do iAngular=1,size(Info_Ang)
    Info_Ang(iAngular)%L_eff = 0
    Info_Ang(iAngular)%nPoints = 0
    if (allocated(Info_Ang(iAngular)%R)) call mma_deallocate(Info_Ang(iAngular)%R)
  end do

end subroutine Close_Info_Ang

subroutine Close_NQ_Data()

  integer iNQ, nNQ

  ! Cleanup and close
  nNQ = size(NQ_data)
  do iNQ=1,nNQ
    call mma_deallocate(NQ_data(iNQ)%Coor)
    NQ_data(iNQ)%A_High = -1.0d99
    NQ_data(iNQ)%A_Low = 1.0d99
    NQ_data(iNQ)%R_RS = 0.0d0
    NQ_data(iNQ)%R_max = 0.0d0
    NQ_data(iNQ)%l_Max = -1
    if (allocated(NQ_data(iNQ)%R_Quad)) call mma_deallocate(NQ_data(iNQ)%R_Quad)
    if (allocated(NQ_data(iNQ)%Angular)) call mma_deallocate(NQ_data(iNQ)%Angular)
    NQ_Data(iNQ)%Atom_Nr = -1
    if (allocated(NQ_data(iNQ)%dOdx)) call mma_deallocate(NQ_data(iNQ)%dOdx)
  end do
  deallocate(NQ_Data)

end subroutine Close_NQ_Data

end module NQ_Structure

