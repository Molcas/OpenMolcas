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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: lRoots, nSPair
real(kind=wp) :: RMat(lRoots**2), ExpX(lRoots**2), X(nSPair)
real(kind=wp), allocatable :: RScr(:)

call mma_allocate(RScr,lRoots**2,Label='RScr')

call ExpMat(ExpX,X,lRoots,nSPair)
call DGEMM_('n','n',lRoots,lRoots,lRoots,One,RMat,lRoots,ExpX,lRoots,Zero,RScr,lRoots)
call DCopy_(lRoots**2,RScr,1,RMat,1)

call mma_deallocate(RScr)

return

end subroutine UpdateRotMat
