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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_XCV_Distrib_SP(mySP,l_mySP,N_mySP)
!
! Thomas Bondo Pedersen, April 2010.
!
! Determine distribution of Shell Pairs according to their
! dimension.

use Para_Info, only: MyRank, nProcs
use Cholesky, only: Cho_Real_Par, nnBstRSh, nnShl, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: l_mySP
integer(kind=iwp), intent(out) :: mySP(l_mySP), N_mySP
integer(kind=iwp) :: iNode, iSP, n
integer(kind=iwp), allocatable :: ProcDim(:)
integer(kind=iwp), external :: Cho_iFindSmallest

#ifdef _DEBUGPRINT_
if (l_mySP < nnShl) call Cho_Quit('Dimension error in Cho_XCV_Distrib_SP',103)
#endif

if (Cho_Real_Par) then
  call mma_allocate(ProcDim,nProcs,Label='ProcDim')
  ProcDim(:) = 0

  N_mySP = 0
  do iSP=1,nnShl
    n = sum(nnBstRSh(1:nSym,iSP,1))
    if (n > 0) then
      iNode = Cho_iFindSmallest(ProcDim,nProcs)
      ProcDim(iNode) = ProcDim(iNode)+n
      if (iNode-1 == myRank) then
        N_mySP = N_mySP+1
        mySP(N_mySP) = iSP
      end if
    end if
  end do

  call mma_deallocate(ProcDim)
else
  N_mySP = 0
  do iSP=1,nnShl
    n = sum(nnBstRSh(1:nSym,iSP,1))
    if (n > 0) then
      N_mySP = N_mySP+1
      mySP(N_mySP) = iSP
    end if
  end do
end if

end subroutine Cho_XCV_Distrib_SP
