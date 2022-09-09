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

module RICD_Info

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: nLen = 10 ! number of elements
integer(kind=iwp) :: iRI_Type = -1
real(kind=wp) :: Thrshld_CD = 1.0e-4_wp
logical(kind=iwp) :: Cho_OneCenter = .false., &
                     Cholesky = .false., &
                     DiagCheck = .false., &
                     Do_acCD_Basis = .true., &
                     Do_RI = .false., &
                     LDF = .false., &
                     LocalDF = .false., &
                     Skip_High_AC = .false.

public :: Cho_OneCenter, Cholesky, DiagCheck, Do_acCD_Basis, Do_RI, iRI_Type, LDF, LocalDF, RICD_Info_Dmp, RICD_Info_Get, &
          Skip_High_AC, Thrshld_CD

contains

subroutine RICD_Info_Dmp()

  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One

  real(kind=wp), allocatable :: rDmp(:)

  call mma_allocate(rDmp,nLen,Label='rDmp:RICD')

  rDmp(01) = real(iRI_Type,kind=wp)
  rDmp(02) = merge(One,Zero,LDF)
  rDmp(03) = merge(One,Zero,Do_RI)
  rDmp(04) = merge(One,Zero,Cholesky)
  rDmp(05) = merge(One,Zero,Do_acCD_Basis)
  rDmp(06) = merge(One,Zero,Skip_High_AC)
  rDmp(07) = merge(One,Zero,Cho_OneCenter)
  rDmp(08) = merge(One,Zero,DiagCheck)
  rDmp(09) = merge(One,Zero,LocalDF)
  rDmp(10) = Thrshld_CD

  call Put_dArray('RICD_Info',rDmp,nLen)
  call mma_deallocate(rDmp)

end subroutine RICD_Info_Dmp

subroutine RICD_Info_Get()

  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero

  real(kind=wp), allocatable :: rDmp(:)

  call mma_allocate(rDmp,nLen,Label='rDmp:RICD')
  call Get_dArray('RICD_Info',rDmp,nLen)

  iRI_Type = nint(rDmp(1))
  LDF = rDmp(2) > Zero
  Do_RI = rDmp(3) > Zero
  Cholesky = rDmp(4) > Zero
  Do_acCD_Basis = rDmp(5) > Zero
  Skip_High_AC = rDmp(6) > Zero
  Cho_OneCenter = rDmp(7) > Zero
  DiagCheck = rDmp(8) > Zero
  LocalDF = rDmp(9) > Zero
  Thrshld_CD = rDmp(10)

  call mma_deallocate(rDmp)

end subroutine RICD_Info_Get

end module RICD_Info
