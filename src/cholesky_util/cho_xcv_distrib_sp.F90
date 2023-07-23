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
use ChoSwp, only: nnBstRSh
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: l_mySP, mySP(l_mySP), N_mySP
#include "cho_para_info.fh"
#include "cholesky.fh"
integer(kind=iwp) :: iNode, iSP, iSym, n
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
    n = nnBstRSh(1,iSP,1)
    do iSym=2,nSym
      n = n+nnBstRSh(iSym,iSP,1)
    end do
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
    n = nnBstRSh(1,iSP,1)
    do iSym=2,nSym
      n = n+nnBstRSh(iSym,iSP,1)
    end do
    if (n > 0) then
      N_mySP = N_mySP+1
      mySP(N_mySP) = iSP
    end if
  end do
end if

end subroutine Cho_XCV_Distrib_SP
