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

subroutine TRPGen(nDim,nAtom,Coor,mTR,CofM,TRVec)

use Slapaf_Info, only: Degen, Smmtrc

implicit real*8(a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
real*8 Coor(3,nAtom), TRVec(3*nAtom*6)
logical CofM
logical, save :: g12K = .true.
real*8, allocatable :: TR(:), Scrt(:), G(:), EVal(:), EVec(:), U(:)

call mma_allocate(TR,18*nAtom,Label='TR')

! Compute the symmetric translations and rotations
!
! B    (nTR x nDim)
!  tr

call TRMake(TR,Coor,nAtom,nTR,Degen,nDim,CofM)

TRVec(1:nTR*nDim) = TR(1:nTR*nDim)

call mma_allocate(Scrt,(3*nAtom)*nTR,Label='Scrt')
call mma_allocate(G,nTR**2,Label='G')
call mma_allocate(EVal,nTR*(nTR+1)/2,Label='EVal')
call mma_allocate(EVec,nTR**2,Label='EVec')

! Eliminate redundancy and produce an orthogonal representation.
!
!    -1/2
! K g        (nTR x mTR)

call mma_allocate(U,nDim,Label='U')
U(:) = One
i = 0
do iAtom=1,nAtom
  do ixyz=1,3
    if (Smmtrc(ixyz,iAtom)) then
      i = i+1
      call DScal_(nTR,sqrt(Degen(ixyz,iAtom)),TRVec((i-1)*nTR+1),1)
    end if
  end do
end do

Thr_ElRed = 1.0D-12
call ElRed(TRVec,nTR,nDim,G,EVal,EVec,mTR,U,Scrt,g12K,Thr_ElRed)

if (mTR > 0) then
  TRVec(1:3*nAtom*nTR) = Zero
  call DGEMM_('T','N',nDim,mTR,nTR,1.0d0,TR,nTR,EVec,nTR,0.0d0,TRVec,nDim)
end if

call mma_deallocate(U)
call mma_deallocate(EVec)
call mma_deallocate(EVal)
call mma_deallocate(G)
call mma_deallocate(Scrt)
call mma_deallocate(TR)

return

end subroutine TRPGen
