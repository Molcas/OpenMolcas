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

use Index_Functions, only: iTri, nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nDim
real(kind=wp), intent(inout) :: V(nDim,nDim)
real(kind=wp), intent(out) :: V12(nDim,nDim)
integer(kind=iwp) :: i, j
real(kind=wp) :: tmp
real(kind=wp), allocatable :: Vec(:,:), VTri(:)

call mma_allocate(Vec,nDim,nDim,Label='Vec')
call mma_allocate(VTri,nTri_Elem(nDim),Label='VTri')

call unitmat(Vec,nDim)

do i=1,nDim
  do j=1,i
    VTri(iTri(i,j)) = V(i,j)
  end do
end do

call JACOB(VTri,Vec,nDim,nDim)

V12(:,:)= Zero
do i=1,nDim
# ifdef _DEBUGPRINT_
  tmp = VTri(iTri(i,i))
  write(u6,*) 'i,tmp=',i,tmp
# endif
  tmp = sqrt(VTri(iTri(i,i)))
  if (tmp < 1.0e-90_wp) then
    V12(i,i) = 1.0e90_wp
  else
    V12(i,i) = One/tmp
  end if
end do

call DGEMM_('N','T',nDim,nDim,nDim,One,V12,nDim,Vec,nDim,Zero,V,nDim)
call DGEMM_('N','N',nDim,nDim,nDim,One,Vec,nDim,V,nDim,Zero,V12,nDim)

call mma_deallocate(VTri)
call mma_deallocate(Vec)

return

end subroutine Compute_V12
