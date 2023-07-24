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

use Index_Functions, only: nTri_Elem
use Slapaf_Info, only: Degen, Smmtrc
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim, nAtom
real(kind=wp), intent(in) :: Coor(3,nAtom)
integer(kind=iwp), intent(out) :: mTR
logical(kind=iwp), intent(in) :: CofM
real(kind=wp), intent(out) :: TRVec(6*3*nAtom)
integer(kind=iwp) :: i, iAtom, ixyz, nTR
real(kind=wp), allocatable :: EVal(:), EVec(:), G(:), Scrt(:), TR(:), U(:)
real(kind=wp), parameter :: Thr_ElRed = 1.0e-12_wp
logical(kind=iwp), parameter :: g12K = .true.

call mma_allocate(TR,6*3*nAtom,Label='TR')

! Compute the symmetric translations and rotations
!
! B    (nTR x nDim)
!  tr

call TRMake(TR,Coor,nAtom,nTR,Degen,nDim,CofM)

TRVec(1:nTR*nDim) = TR(1:nTR*nDim)

call mma_allocate(Scrt,(3*nAtom)*nTR,Label='Scrt')
call mma_allocate(G,nTR**2,Label='G')
call mma_allocate(EVal,nTri_Elem(nTR),Label='EVal')
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
      TRVec((i-1)*nTR+1:i*nTR) = sqrt(Degen(ixyz,iAtom))*TRVec((i-1)*nTR+1:i*nTR)
    end if
  end do
end do

call ElRed(TRVec,nTR,nDim,G,EVal,EVec,mTR,U,Scrt,g12K,Thr_ElRed)

if (mTR > 0) call DGEMM_('T','N',nDim,mTR,nTR,One,TR,nTR,EVec,nTR,Zero,TRVec,nDim)

call mma_deallocate(U)
call mma_deallocate(EVec)
call mma_deallocate(EVal)
call mma_deallocate(G)
call mma_deallocate(Scrt)
call mma_deallocate(TR)

return

end subroutine TRPGen
