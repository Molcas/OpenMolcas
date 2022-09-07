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

subroutine RI_XDiag(Diag,nDiag)
! Thomas Bondo Pedersen, Jan. 2007.
!
! Purpose: compute exact integral diagonal.

use ChoArr, only: iSP2F, nBstSh
use ChoSwp, only: iiBstRSh, IndRed, nnBstRSh
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nDiag
real(kind=wp) :: Diag(nDiag)
#include "cholesky.fh"
integer(kind=iwp) :: i, i1, i2, ID, iSAB, iShlA, iShlB, iSym, l_SewMem, NumAB
real(kind=wp), allocatable :: Scr(:)
logical(kind=iwp), external :: Rsv_Tsk

! Allocate memory.
! ----------------

call Init_Tsk(ID,nnShl)

call mma_allocate(Scr,Mx2Sh,Label='Scr')
call mma_maxDBLE(l_SewMem)

! Initialize diagonal array.
! --------------------------

call fZero(Diag,nnBstRT(1))

! Parallel loop over shell pairs in first red. set.
! -------------------------------------------------

do while (Rsv_Tsk(ID,iSAB))

  ! Get shells.
  ! -----------

  call Cho_InvPck(iSP2F(iSAB),iShlA,iShlB,.true.)

  ! Compute (AB|AB).
  ! ----------------

  if (iShlA == iShlB) then
    NumAB = nBstSh(iShlA)*(nBstSh(iShlA)+1)/2
  else
    NumAB = nBstSh(iShlA)*nBstSh(iShlB)
  end if
  ShA = iShlA
  ShB = iShlB
  call Cho_MCA_DiagInt(iShlA,iShlB,Scr,NumAB)

  ! Extract diagonal elements.
  ! --------------------------

  do iSym=1,nSym
    i1 = iiBstR(iSym,1)+iiBstRSh(iSym,iSAB,1)+1
    i2 = i1+nnBstRSh(iSym,iSAB,1)-1
    do i=i1,i2
      Diag(i) = Scr(IndRed(i,1))
    end do
  end do

end do
call Free_Tsk(ID)

! Sync diagonal.
! --------------

call GAdGOP(Diag,nnBstRT(1),'+')

! Deallocate memory.
! ------------------

call xRlsMem_Ints()
call mma_deallocate(Scr)

end subroutine RI_XDiag
