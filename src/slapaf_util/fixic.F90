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

subroutine FIXIC(nFix,SS,mInt,B,NDIM,F,Label,u)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nFix, mInt, NDIM
real(kind=wp), intent(inout) :: SS(mInt)
real(kind=wp), intent(in) :: B(nDim*mInt), u(nDim)
real(kind=wp), intent(out) :: F(nDim)
character(len=8), intent(in) :: Label(mInt)
integer(kind=iwp) :: I
real(kind=wp), allocatable :: uB(:,:), uInv(:,:)

! write out the internal coordinates which will be fixed

write(u6,*)
write(u6,*) ' Following internal coordinates are fixed'
write(u6,*)

! loop over all internal coordinates to be fixed

do I=mInt-nFix+1,mInt
  write(u6,'(A,A,ES10.3,A)') Label(i),' with a gradient of ',SS(I),' is frozen and the gradient is annihilated'
  SS(i) = Zero
end do

! now transform remaining internal coordinates back to cartesian ba
!                      -1 +
!                fx = u  B  fq

call mma_allocate(uInv,nDim,nDim,Label='uInv')
uInv(:,:) = Zero
do i=1,nDim
  uInv(i,i) = One/u(i)
end do
!call RecPrt('uInv',' ',uInv,nDim,nDim)

call mma_allocate(uB,mInt,nDim,Label='uB')
uB(:,:) = Zero

call DGEMM_('N','N',nDim,mInt,nDim,One,uInv,nDim,B,nDim,Zero,uB,nDim)
!call RecPrt('uInvB',' ',uB,nDim,mInt)
call DGEMM_('N','N',nDim,1,mInt,One,uB,nDim,SS,mInt,Zero,F,nDim)
!call RecPrt('F',' ',F,mInt,1)

call mma_deallocate(uB)
call mma_deallocate(uInv)

return

end subroutine FIXIC
