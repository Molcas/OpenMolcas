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
! Copyright (C) 2017, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Remove_TR(nQ,nX,nQQ,KMat,nK,TRVec,nTR,BM,iBM,nqB,nB)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nQ, nX, nQQ, nK, nTR, nB, iBM(nB), nqB(nQ)
real(kind=wp) :: KMat(nQ,nK), TRVec(nX,nTR), BM(nB)
integer(kind=iwp) :: i, iB, iK, iQ, iTR, iV, iX
real(kind=wp), allocatable :: TR(:), VecInt(:)
real(kind=wp), external :: DDot_

call mma_allocate(TR,nK)
call mma_allocate(VecInt,nX)

call FZero(TR,nK)

do iK=1,nK

  ! Get the normalized Cartesian vector for this internal coordinate
  call FZero(VecInt,nX)
  iB = 0
  do iQ=1,nQ
    do i=1,nqB(iQ)
      iB = iB+1
      iX = iBM(iB)
      VecInt(iX) = VecInt(iX)+KMat(iQ,iK)*BM(iB)
    end do
  end do
  call dScal_(nX,One/sqrt(dDot_(nX,VecInt,1,VecInt,1)),VecInt,1)

  ! Compute the overlap with the external translations and rotations
  do iTR=1,nTR
    TR(iK) = TR(iK)+dDot_(nX,VecInt,1,TRvec(1,iTR),1)**2
  end do
end do

! Put the nK-nQQ vectors with largest overlap at the end
do iK=nK,nQQ+1,-1
  iV = maxloc(TR(1:iK),1)
  if (iV /= iK) call dSwap_(nQ,KMat(1,iK),1,KMat(1,iV),1)
end do

call mma_deallocate(TR)
call mma_deallocate(VecInt)

return

end subroutine Remove_TR
