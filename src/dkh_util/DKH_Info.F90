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

module DKH_Info

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: nCtrLD = 0, iCtrLD(10) = [0,0,0,0,0,0,0,0,0,0], iRelae_Info = 0
real(kind=wp) :: radiLD = Zero
logical(kind=iwp) :: DKroll = .false., LDKroll = .false., BSS = .false.

public :: nCtrLD, iCtrLD, radiLD, DKroll, LDKroll, BSS, DKH_Info_Get, DKH_Info_Dmp

contains

subroutine DKH_Info_Dmp()
  use stdalloc, only: mma_allocate, mma_deallocate
  real(kind=wp), allocatable :: rDmp(:)
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: Length = 1+10+5
# include "relae.fh"

  iRelae_Info = iRelae

  call mma_allocate(rDmp,Length,Label='rDmp:DKH')
  rDmp(1) = real(nCtrLD,kind=wp)
  do i=1,10
    rDmp(i+1) = iCtrLD(i)
  end do
  rDmp(12) = radiLD
  rDmp(13) = merge(One,Zero,DKroll)
  rDmp(14) = merge(One,Zero,LDKroll)
  rDmp(15) = merge(One,Zero,BSS)
  rDmp(16) = real(iRelae_Info,kind=wp)
  call Put_dArray('DKH_Info',rDmp,Length)
  call mma_deallocate(rDmp)
end subroutine DKH_Info_Dmp

subroutine DKH_Info_Get()
  use stdalloc, only: mma_allocate, mma_deallocate
  real(kind=wp), allocatable :: rDmp(:)
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: Length = 1+10+5
# include "relae.fh"

  call mma_allocate(rDmp,Length,Label='rDmp:DKH')
  call Get_dArray('DKH_Info',rDmp,Length)

  nCtrLD = nint(rDmp(1))
  do i=1,10
    iCtrLD(i) = nint(rDmp(i+1))
  end do
  radiLD = rDmp(12)
  DKroll = nint(rDmp(13)) == 1
  LDKroll = nint(rDmp(14)) == 1
  BSS = nint(rDmp(15)) == 1
  iRelae_Info = nint(rDmp(16))

  iRelae = iRelae_Info

  call mma_deallocate(rDmp)
end subroutine DKH_Info_Get

end module DKH_Info
