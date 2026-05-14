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
! Copyright (C) 2007, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_P_Distrib_SP_byDim(mySP,N_mySP)
!
! Thomas Bondo Pedersen, June 2007.
!
! Determine distribution of ShellPairs by dimension.

use Index_Functions, only: nTri_Elem
use Para_Info, only: MyRank, nProcs
use Cholesky, only: Cho_Real_Par, iSP2F, nBstSh, nnShl
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(_OUT_) :: mySP(*)
integer(kind=iwp), intent(out) :: N_mySP
integer(kind=iwp) :: iA, iAB, iB, iNode, iSP, n
integer(kind=iwp), allocatable :: ProcDim(:)
integer(kind=iwp), external :: Cho_iFindSmallest

if (Cho_Real_Par) then

  call mma_allocate(ProcDim,[0,nProcs-1],Label='ProcDim')
  ProcDim(:) = 0

  N_mySP = 0
  do iSP=1,nnShl
    iAB = iSP2F(iSP)
    call Cho_InvPck(iAB,iA,iB,.true.)
    if (iA == iB) then
      n = nTri_Elem(nBstSh(iA))
    else
      n = nBstSh(iA)*nBstSh(iB)
    end if
    iNode = Cho_iFindSmallest(ProcDim,size(ProcDim))-1
    ProcDim(iNode) = ProcDim(iNode)+n
    if (iNode == myRank) then
      N_mySP = N_mySP+1
      mySP(N_mySP) = iSP
    end if
  end do

  call mma_deallocate(ProcDim)

else

  N_mySP = nnShl
  do iSP=1,N_mySP
    mySP(iSP) = iSP
  end do

end if

end subroutine Cho_P_Distrib_SP_byDim
