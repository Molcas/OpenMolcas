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

subroutine Reset_ThrGrd(nIter,mTtAtm,ThrGrd)

use Slapaf_Info, only: Cx, nDimBC
use Slapaf_procedures, only: Box, Hidden
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Ten, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nIter, mTtAtm
real(kind=wp), intent(inout) :: ThrGrd
integer(kind=iwp) :: i, iIter, mTR, nBonds, nHidden, nMax, nSaddle, nsAtom
logical(kind=iwp) :: Found
integer(kind=iwp), allocatable :: AN(:), TabA(:,:,:), TabAI(:), TabB(:,:)
real(kind=wp), allocatable :: Coor(:,:), Tmp(:), TR(:), Vec(:)

!                                                                      *
!***********************************************************************
!                                                                      *
nsAtom = size(Cx,2)
!                                                                      *
!***********************************************************************
!                                                                      *
call qpg_dArray('Saddle',Found,nSaddle)
if (.not. Found) return
iIter = nIter           ! Normal Computation
!                                                                      *
!***********************************************************************
!                                                                      *
! Find the translational and rotational eigenvectors for the
! current structure.

call mma_allocate(TR,18*nsAtom,Label='TR')
TR(:) = Zero

call TRPGen(nDimBC,nsAtom,Cx(1,1,iIter),mTR,.false.,TR)

!call RecPrt('TR',' ',TR,nDimBC,mTR)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_Allocate(TabAI,2*mTtAtm,Label='TabAI')
call mma_Allocate(Vec,3*mTtAtm*nDimBC,Label='Vec')
call mma_Allocate(AN,mTtAtm,Label='AN')
call mma_Allocate(Coor,3,mTtAtm,Label='Coor')

! Generate Grand atoms list

call GenCoo(Cx(1,1,iIter),nsAtom,Coor,mTtAtm,Vec,nDimBC,AN,TabAI)

! Are there some hidden frozen atoms ?

call Hidden(Coor,AN,nHidden)

! Generate bond list

call Box(Coor,mTtAtm+nHidden,AN,TabB,TabA,nBonds,nMax)
!                                                                      *
!***********************************************************************
!                                                                      *
! If there are some bond types 2, and we are in saddle
! far from the TS, let us get a reduced threshold to avoid
! wasting our time

call mma_allocate(Tmp,nSaddle,Label='Tmp')
call Get_dArray('Saddle',Tmp,nSaddle)
Found = .false.
if (Tmp(nSaddle-1) > Half) then
  do i=1,nBonds
    if (TabB(3,i) == 2) then
      Found = .true.
      exit
    end if
  end do
  if (Found) then
    !ThrGrd = 0.03_wp
    ThrGrd = Ten*ThrGrd
    call WarningMessage(1,'Molecule composed of many fragments Convergence threshold reduced')
  end if
end if
call mma_deallocate(Tmp)

call mma_deallocate(TabA)
call mma_deallocate(TabB)
call mma_deallocate(Coor)
call mma_deallocate(AN)
call mma_deallocate(Vec)
call mma_deallocate(TabAI)
call mma_deallocate(TR)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Reset_ThrGrd
