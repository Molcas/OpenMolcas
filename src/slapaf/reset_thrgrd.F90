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

use Slapaf_Info, only: Cx
use Slapaf_Parameters, only: nDimBC

implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
logical Found
integer, allocatable :: TabAI(:), AN(:)
real*8, allocatable :: TR(:), Vec(:), Coor(:,:), Tmp(:)
integer, allocatable :: TabB(:,:), TabA(:,:,:)
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine Box(Coor,nsAtom,iANr,TabB,TabA,nBonds,nMax)
    integer nsAtom
    real*8 Coor(3,nsAtom)
    integer iANr(nsAtom)
    integer, allocatable :: TabB(:,:), TabA(:,:,:)
    integer nBonds, nMax
  end subroutine Box
  subroutine Hidden(Coor,AN,nHidden)
    real*8, allocatable :: Coor(:,:)
    integer, allocatable :: AN(:)
    integer nHidden
  end subroutine Hidden
end interface

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

mTtAtm = mTtAtm+nHidden
call Box(Coor,mTtAtm,AN,TabB,TabA,nBonds,nMax)
mTtAtm = mTtAtm-nHidden
!                                                                      *
!***********************************************************************
!                                                                      *
! If there are some bond types 2, and we are in saddle
! far from the TS, let us get a reduced threshold to avoid
! wasting our time

call mma_allocate(Tmp,nSaddle,Label='Tmp')
call Get_dArray('Saddle',Tmp,nSaddle)
Found = .false.
if (Tmp(nSaddle-1) > 0.50d0) then
  do i=1,nBonds
    if (TabB(3,i) == 2) then
      Found = .true.
      Go To 20
    end if
  end do
20 continue
  if (Found) then
    ! ThrGrd = 0.03D0
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
