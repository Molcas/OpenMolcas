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

subroutine Cho_VecDsk_GetLQ(QVec,l_QVec,LstQSP,nQSP,iV1,nV,mSym)
!
! Purpose: extract elements corresponding to qualified columns of
!          vectors on disk.

use Cholesky, only: nnBstRSh, nQual, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: l_QVec, nQSP, LstQSP(nQSP), mSym, iV1(mSym), nV(mSym)
real(kind=wp), intent(inout) :: QVec(l_QVec)
integer(kind=iwp) :: iAB, iLoc, iQ, iQSP, iRedC, iRedQ, iShlAB, iSym, iV2, jLoc, jV, kOffQ, kQ, l_iQuAB_2, l_Scr, lDim, &
                     nVecTot(8), nVT
integer(kind=iwp), allocatable :: iQuAB_2(:)
real(kind=wp), allocatable :: Scr(:)
character(len=*), parameter :: SecNam = 'Cho_VecDsk_GetLQ'
integer(kind=iwp), external :: Cho_P_LocalSP

! Check input.
! ------------

if (nQSP < 1) return

if (mSym < nSym) call Cho_Quit('mSym<nSym in '//SecNam,104)

nVT = sum(nV(1:nSym))
if (nVT < 1) return

do iSym=1,nSym
  if (nV(iSym) > 0) then
    if (iV1(iSym) < 1) call Cho_Quit('iV1<1 in '//SecNam,104)
  end if
end do

! Allocate memory for reading.
! ----------------------------

l_Scr = 0
do iSym=1,nSym
  if ((nV(iSym) > 0) .and. (nQual(iSym) > 0)) then
    lDim = 0
    do iQSP=1,nQSP
      iShlAB = Cho_P_LocalSP(LstQSP(iQSP))
      lDim = lDim+nnBstRSh(iSym,iShlAB,1)
    end do
    l_Scr = max(l_Scr,lDim)
  end if
end do
call mma_allocate(Scr,l_Scr,Label='Scr')

! Allocate extra mapping from qualified to reduced set.
! -----------------------------------------------------

l_iQuAB_2 = nQual(1)
do iSym=2,nSym
  l_iQuAB_2 = max(l_iQuAB_2,nQual(iSym))
end do
call mma_allocate(iQuAB_2,l_iQuAB_2,Label='iQuAB_2')

! Extract in each symmetry block.
! -------------------------------

call Cho_P_GetGV(nVecTot,nSym)

jLoc = 2
iLoc = 3
iRedC = -1
kOffQ = 0
do iSym=1,nSym
  iRedQ = -2
  if ((nV(iSym) > 0) .and. (nQual(iSym) > 0)) then
    iV2 = iV1(iSym)+nV(iSym)-1
    do jV=iV1(iSym),iV2
      call Cho_1VecRd_SP(Scr,l_Scr,jV,iSym,LstQSP,nQSP,iRedC,iLoc)
      if (iRedQ /= iRedC) then
        call Cho_SetQ2(iQuAB_2,LstQSP,nQSP,iSym,jLoc,iLoc)
        iRedQ = iRedC
      end if
      kQ = kOffQ+nQual(iSym)*(jV-1)
      do iQ=1,nQual(iSym)
        iAB = iQuAB_2(iQ)
        QVec(kQ+iQ) = Scr(iAB)
      end do
    end do
  end if
  kOffQ = kOffQ+nQual(iSym)*nVecTot(iSym)
end do

! Deallocations.
! --------------

call mma_deallocate(iQuAB_2)
call mma_deallocate(Scr)

end subroutine Cho_VecDsk_GetLQ
