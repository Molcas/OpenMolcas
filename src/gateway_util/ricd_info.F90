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

private
public :: iRI_Type, LDF, Do_RI, Cholesky, Do_acCD_Basis, Skip_High_AC, &
          Cho_OneCenter, DiagCheck, LocalDF, Do_nacCD_Basis, Thrshld_CD, &
          RICD_Info_Dmp, RICD_Info_Get

integer :: iRI_Type = -1
logical :: LDF = .false.
logical :: Do_RI = .false.
logical :: Cholesky = .false.
logical :: Do_acCD_Basis = .true.
logical :: Skip_High_AC = .false.
logical :: Cho_OneCenter = .false.
logical :: DiagCheck = .false.
logical :: LocalDF = .false.
logical :: Do_nacCD_Basis = .false.
real*8 :: Thrshld_CD = 1.0D-4
#include "stdalloc.fh"

contains

subroutine RICD_Info_Dmp()

  real*8, allocatable :: rDmp(:)
  integer i
  integer :: Len = 11

  call mma_allocate(rDmp,Len,Label='rDmp:RICD')

  rDmp(1) = dble(iRI_Type)
  i = 0
  if (LDF) i = 1
  rDmp(2) = dble(i)
  i = 0
  if (Do_RI) i = 1
  rDmp(3) = dble(i)
  i = 0
  if (Cholesky) i = 1
  rDmp(4) = dble(i)
  i = 0
  if (Do_acCD_Basis) i = 1
  rDmp(5) = dble(i)
  i = 0
  if (Skip_High_AC) i = 1
  rDmp(6) = dble(i)
  i = 0
  if (Cho_OneCenter) i = 1
  rDmp(7) = dble(i)
  i = 0
  if (DiagCheck) i = 1
  rDmp(8) = dble(i)
  i = 0
  if (LocalDF) i = 1
  rDmp(9) = dble(i)
  i = 0
  if (Do_nacCD_Basis) i = 1
  rDmp(10) = dble(i)
  rDmp(11) = Thrshld_CD

  call Put_dArray('RICD_Info',rDmp,Len)
  call mma_deallocate(rDmp)

end subroutine RICD_Info_Dmp

subroutine RICD_Info_Get()

  real*8, allocatable :: rDmp(:)
  integer :: Len = 11

  call mma_allocate(rDmp,Len,Label='rDmp:RICD')
  call Get_dArray('RICD_Info',rDmp,Len)

  iRI_Type = nint(rDmp(1))
  LDF = nint(rDmp(2)) == 1
  Do_RI = nint(rDmp(3)) == 1
  Cholesky = nint(rDmp(4)) == 1
  Do_acCD_Basis = nint(rDmp(5)) == 1
  Skip_High_AC = nint(rDmp(6)) == 1
  Cho_OneCenter = nint(rDmp(7)) == 1
  DiagCheck = nint(rDmp(8)) == 1
  LocalDF = nint(rDmp(9)) == 1
  Do_nacCD_Basis = nint(rDmp(10)) == 1
  Thrshld_CD = rDmp(11)

  call mma_deallocate(rDmp)

end subroutine RICD_Info_Get

end module RICD_Info
