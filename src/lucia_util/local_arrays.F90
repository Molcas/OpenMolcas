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

module Local_Arrays

! CLBT  : Length of each Batch (in blocks)
! CLEBT : Length of each Batch (in elements)
! CI1BT : Length of each block
! CIBT  : Info on each block
! CBLTP : BLock type for each symmetry

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
private

integer(kind=iwp), allocatable :: CI1BT(:), CIBT(:), CLBT(:), CLEBT(:)
integer(kind=iwp), allocatable, target :: CBLTP(:)

public :: Allocate_Local_Arrays, CBLTP, CI1BT, CIBT, CLBT, CLEBT, Deallocate_Local_Arrays

contains

subroutine Allocate_Local_Arrays(MXNTTS,NSMST)

  integer(kind=iwp), intent(in) :: MXNTTS, NSMST

  if (allocated(CLBT)) call Abend()
  call mma_allocate(CLBT,MXNTTS,Label='CLBT',safe='*')
  call mma_allocate(CLEBT,MXNTTS,Label='CLEBT',safe='*')
  call mma_allocate(CI1BT,MXNTTS,Label='CI1BT',safe='*')
  call mma_allocate(CIBT,8*MXNTTS,Label='CIBT',safe='*')
  call mma_allocate(CBLTP,NSMST,Label='CBLTP',safe='*')

end subroutine Allocate_Local_Arrays

subroutine Deallocate_Local_Arrays()

  call mma_deallocate(CLBT,safe='*')
  call mma_deallocate(CLEBT,safe='*')
  call mma_deallocate(CI1BT,safe='*')
  call mma_deallocate(CIBT,safe='*')
  call mma_deallocate(CBLTP,safe='*')

end subroutine Deallocate_Local_Arrays

end module Local_Arrays
