!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 07, 2022, created this file.               *
!*****************************************************************

subroutine UpdateRotMat(RMat,ExpX,X,lRoots,nSPair)

use Index_Functions, only: iTri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lRoots, nSPair
real(kind=wp), intent(inout) :: RMat(lRoots**2)
real(kind=wp), intent(out) :: ExpX(lRoots,lRoots)
real(kind=wp), intent(in) :: X(nSPair)
integer(kind=iwp) :: I, iIJ, J
real(kind=wp) :: maxtheta
real(kind=wp), allocatable :: RScr(:)

call mma_allocate(RScr,lRoots**2,Label='RScr')

ExpX(1,1) = Zero
do I=2,lRoots
  ExpX(I,I) = Zero
  do J=1,I-1
    iIJ = iTri(I-1,J)
    ExpX(J,I) = X(iIJ)
    ExpX(I,J) = -X(iIJ)
  end do
end do
call Exp_eig(lRoots,ExpX,maxtheta)
call DGEMM_('n','n',lRoots,lRoots,lRoots,One,RMat,lRoots,ExpX,lRoots,Zero,RScr,lRoots)
RMat(:) = RScr(:)

call mma_deallocate(RScr)

return

end subroutine UpdateRotMat
