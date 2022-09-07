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

subroutine Compute_V12(V,V12,nDim)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nDim
real(kind=wp) :: V(nDim,nDim), V12(nDim,nDim)
real(kind=wp), allocatable :: Vec(:), VTri(:)

call mma_allocate(Vec,nDim**2,Label='Vec')
call mma_allocate(VTri,nDim*(nDim+1)/2,Label='VTri')

call Compute_V12_Inner(V,V12,VTri,Vec,nDim)

call mma_deallocate(VTri)
call mma_deallocate(Vec)

return

end subroutine Compute_V12
