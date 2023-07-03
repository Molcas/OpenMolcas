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

implicit real*8(a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
real*8 SS(mInt), B(nDim*mInt), F(nDim), u(nDim)
character(len=8) Label(mInt)
real*8, allocatable :: uInv(:,:), uB(:,:)

! write out the internal coordinates which will be fixed

write(6,*)
write(6,*) ' Following internal coordinates are fixed'
write(6,*)

! loop over all internal coordinates to be fixed

do I=mInt-nFix+1,mInt
  write(6,'(A,A,E10.3,A)') Label(i),' with a gradient of ',SS(I),' is frozen and the gradient is annihilated'
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
