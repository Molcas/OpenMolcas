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

use Constants, only: Zero, One, c_in_au
use Definitions, only: wp, iwp

implicit none
private

! IRELMP = 0  .... NOPAIR (DK2)
!        = 1  .... NOPAIR (DK1)
!        = 2  .... NOPAIR (DK2)
!        = 3  .... NOPAIR (DK3)
!        = 4  .... full NOPAIR (DK3)
!        = 11 .... RESC
!        = 21 .... ZORA
!        = 22 .... ZORA-FP
!        = 23 .... IORA
!
! IRELAE = 0  .... DKH
!        = 1  .... DK1
!        = 2  .... DK2
!        = 3  .... DK3
!        = 4  .... DK3full
!        = 11 .... RESC
!        = 21 .... ZORA
!        = 22 .... ZORA(FP)
!        = 23 .... IORA
!        = 101.... X2C
!        = 102.... BSS
!
! NB: The IRELAE flag has been extended to account for
!     arbitrary-order DKH with different parametrizations!
!     IMPORTANT: new arbitrary-order DKH routines are only
!                called for IRELAE values LARGER than 1000.

integer(kind=iwp) :: nCtrLD = 0, iCtrLD(10) = [0,0,0,0,0,0,0,0,0,0], iRelae = -1, iRelmp = -1
real(kind=wp) :: radiLD = Zero, cLightAU = c_in_au
logical(kind=iwp) :: DKroll = .false., LDKroll = .false., BSS = .false.

public :: nCtrLD, iCtrLD, radiLD, DKroll, LDKroll, BSS, iRelae, iRelmp, cLightAU, DKH_Info_Get, DKH_Info_Dmp

contains

subroutine DKH_Info_Dmp()
  use stdalloc, only: mma_allocate, mma_deallocate
  real(kind=wp), allocatable :: rDmp(:)
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: Length = 17

  call mma_allocate(rDmp,Length,Label='rDmp:DKH')
  rDmp(1) = real(nCtrLD,kind=wp)
  do i=1,10
    rDmp(i+1) = iCtrLD(i)
  end do
  rDmp(12) = radiLD
  rDmp(13) = merge(One,Zero,DKroll)
  rDmp(14) = merge(One,Zero,LDKroll)
  rDmp(15) = merge(One,Zero,BSS)
  rDmp(16) = cLightAU
  rDmp(17) = real(iRelae,kind=wp)
  call Put_dArray('DKH_Info',rDmp,Length)
  call mma_deallocate(rDmp)
end subroutine DKH_Info_Dmp

subroutine DKH_Info_Get()
  use stdalloc, only: mma_allocate, mma_deallocate
  real(kind=wp), allocatable :: rDmp(:)
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: Length = 17

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
  cLightAU = rDmp(16)
  iRelae = nint(rDmp(17))

  call mma_deallocate(rDmp)
end subroutine DKH_Info_Get

end module DKH_Info
