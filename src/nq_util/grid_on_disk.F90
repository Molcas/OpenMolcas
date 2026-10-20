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

module Grid_On_Disk

use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: Use_Old = 0, Regenerate = 1, &
                                Not_Specified = 0, Final_Grid = 1, Intermediate = 2
integer(kind=iwp) :: G_S(2), Grid_Status, iDisk_Grid, iDisk_Set(2), iGrid_Set, Lu_Grid, LuGridFile, nBatch, nBatch_Max = 128, &
                     Old_Functional_Type
integer(kind=iwp), allocatable :: GridInfo(:,:), iBatchInfo(:,:)

public :: ExpandBatchInfo, Final_Grid, G_S, Grid_Status, GridInfo, iBatchInfo, iDisk_Grid, iDisk_Set, iGrid_Set, Intermediate, &
          Lu_Grid, LuGridFile, nBatch, nBatch_Max, Not_Specified, Old_Functional_Type, Regenerate, Use_Old

contains

subroutine ExpandBatchInfo()

  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp) :: n
  integer(kind=iwp), allocatable :: new_iBatchInfo(:,:)

  n = size(iBatchInfo,2)
  nBatch_Max = 2*n
  call mma_allocate(new_iBatchInfo,size(iBatchInfo,1),nBatch_Max,label='new_iBatchInfo')
  new_iBatchInfo(:,1:n) = iBatchInfo
  new_iBatchInfo(:,n+1:) = 0
  call mma_deallocate(iBatchInfo)
  call move_alloc(new_iBatchInfo,iBatchInfo)

end subroutine ExpandBatchInfo

end module Grid_On_Disk
