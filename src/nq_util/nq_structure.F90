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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: LMax_NQ = 62

type Info_Ang_t
  integer(kind=iwp) :: L_eff = 0
  integer(kind=iwp) :: nPoints = 0
  real(kind=wp), allocatable :: R(:,:)
end type Info_Ang_t

type NQ_data_t
  real(kind=wp) :: Coor(3) = Zero
  real(kind=wp) :: A_High = -huge(Zero)
  real(kind=wp) :: A_Low = huge(Zero)
  real(kind=wp) :: R_RS = Zero
  real(kind=wp) :: R_max = Zero
  integer(kind=iwp) :: l_max = -1
  real(kind=wp), allocatable :: R_Quad(:,:)
  integer(kind=iwp), allocatable :: Angular(:)
  integer(kind=iwp) :: Atom_Nr = -1
  real(kind=wp) :: dOdx(3,3,3) = Zero
  integer(kind=iwp) :: Grad_idx(3) = 0
  integer(kind=iwp) :: Fact = 0
end type NQ_data_t

type(Info_Ang_t) Info_Ang(LMax_NQ)
type(NQ_data_t), allocatable :: NQ_data(:)

public :: Close_Info_Ang, Close_NQ_Data, Info_Ang, LMax_NQ, NQ_data, Open_NQ_Data

! Private extensions to mma interfaces

interface mma_allocate
  module procedure :: nqdata_mma_allo_1D, nqdata_mma_allo_1D_lim
end interface
interface mma_deallocate
  module procedure :: nqdata_mma_free_1D
end interface

contains

subroutine Open_NQ_Data(Coor)

  real(kind=wp), intent(in) :: Coor(:,:)
  integer(kind=iwp) :: iNQ, nNQ

  nNQ = size(Coor,2)
  call mma_allocate(NQ_data,nNQ,'NQ_data')
  do iNQ=1,nNQ
    NQ_data(iNQ)%Coor(:) = Coor(1:3,iNQ)
  end do

# include "macros.fh"
  unused_proc(mma_allocate(NQ_data,[0,0]))

end subroutine Open_NQ_Data

subroutine Close_Info_Ang()

  integer(kind=iwp) :: iAngular

  do iAngular=1,size(Info_Ang)
    Info_Ang(iAngular)%L_eff = 0
    Info_Ang(iAngular)%nPoints = 0
    call mma_deallocate(Info_Ang(iAngular)%R,safe='*')
  end do

end subroutine Close_Info_Ang

subroutine Close_NQ_Data()

  integer(kind=iwp) :: iNQ

  ! Cleanup and close
  do iNQ=1,size(NQ_data)
    call mma_deallocate(NQ_data(iNQ)%R_Quad,safe='*')
    call mma_deallocate(NQ_data(iNQ)%Angular,safe='*')
  end do
  call mma_deallocate(NQ_Data)

end subroutine Close_NQ_Data

! Private extensions to mma_interfaces, using preprocessor templates
! (see src/mma_util/stdalloc.f)

! Define nqdata_mma_allo_1D, nqdata_mma_allo_1D_lim, nqdata_mma_free_1D
! (using _NO_GARBLE_ because all members are initialized)
#define _TYPE_ type(NQ_data_t)
#  define _NO_GARBLE_
#  define _SUBR_NAME_ nqdata_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'nqd_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

end module NQ_Structure
