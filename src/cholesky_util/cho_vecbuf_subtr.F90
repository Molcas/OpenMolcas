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

subroutine Cho_VecBuf_Subtr(xInt,Wrk,lWrk,iSym,DoTime,DoStat)
!
! Purpose: subtract contributions to qualified columns from the
!          vectors stored in the buffer (if any).
!
! DoTime: time as vector subtraction.
! DpStat: update statistics info (#calls to dGeMM).

use ChoSwp, only: iQuAB, nnBstRSh, iiBstRSh
use ChoArr, only: LQ
use ChoVecBuf, only: CHVBUF, ip_CHVBUF_SYM, l_CHVBUF_SYM, nVec_in_Buf
use ChoSubScr, only: Cho_SScreen, SSTau, SubScrStat, DSubScr, DSPNm, SSNorm

implicit real*8(a-h,o-z)
real*8, target :: xInt(*), Wrk(lWrk)
logical DoTime, DoStat
#include "cholesky.fh"
character(len=16), parameter :: SecNam = 'Cho_VecBuf_Subtr'
#ifdef _DEBUGPRINT_
logical, parameter :: LocDbg = .true.
#else
logical, parameter :: LocDbg = .false.
#endif
real*8, parameter :: xMOne = -1.0d0, One = 1.0d0
real*8, pointer :: V(:,:) => null(), U(:,:) => null(), W(:,:) => null()

! Return if nothing to do.
! ------------------------

if (l_ChVBuf_Sym(iSym) < 1) then
  if (LocDbg) then
    write(Lupri,*) SecNam,': returns immediately!'
    write(Lupri,*) ' -- no buffer allocated for sym. ',iSym
  end if
  return
end if
if (nVec_in_Buf(iSym) < 1) then
  if (LocDbg) then
    write(Lupri,*) SecNam,': returns immediately!'
    write(Lupri,*) ' -- buffer is empty for sym. ',iSym
  end if
  return
end if
if (nQual(iSym) < 1) then
  if (LocDbg) then
    write(Lupri,*) SecNam,': returns immediately!'
    write(Lupri,*) ' -- no qualified columns of sym. ',iSym
  end if
  return
end if
if (nnBstR(iSym,2) < 1) then
  if (LocDbg) then
    write(Lupri,*) SecNam,': returns immediately!'
    write(Lupri,*) ' -- empty symmetry block (sym. ',iSym,')'
  end if
  return
end if

! Start timing.
! -------------

if (DoTime) call Cho_Timer(C1,W1)

! Initialize.
! -----------

xTot = 0.0d0
xDon = 0.0d0

! Set up vector batch.
! --------------------

nVec = min(lWrk/nQual(iSym),nVec_in_Buf(iSym))
if (nVec < 1) then
  call Cho_Quit('Insufficient memory for batch in '//SecNam,101)
  nBatch = -999999 ! avoid compiler warnings
else
  nBatch = (nVec_in_Buf(iSym)-1)/nVec+1
end if

! Map the integral array, xInt, onto the pointer U
! ------------------------------------------------

lRow = nnBstR(iSym,2)
lCol = nQual(iSym)
iS = 1
iE = iS-1+lRow*lCol

U(1:lRow,1:lCol) => xInt(iS:iE)

! Start batch loop.
! -----------------

do iBatch=1,nBatch

  ! Set info for this batch.
  ! ------------------------

  if (iBatch == nBatch) then
    NumV = nVec_in_Buf(iSym)-nVec*(nBatch-1)
  else
    NumV = nVec
  end if
  iVec0 = nVec*(iBatch-1)

# ifdef _DEBUGPRINT_
  Need = nQual(iSym)*NumV
  if (lWrk < Need) call Cho_Quit('Batch setup error in '//SecNam,104)
# endif

  lRow = nnBstR(iSym,2)
  lCol = iVec0+NumV
  iS = ip_ChVBuf_Sym(iSym)
  iE = iS-1+lRow*lCol

  V(1:lRow,1:lCol) => CHVBUF(iS:iE)

  ! Screened or unscreened subtraction section.
  ! The screened version uses level 2 blas, while the unscreened
  ! one employs level 3 blas.
  ! ------------------------------------------------------------

  if (Cho_SScreen) then

    lRow = NumV
    lCol = nQual(iSym)
    iS = 1
    iE = iS-1+lRow*lCol

    W(1:lRow,1:lCol) => Wrk(iS:iE)

    ! Copy out sub-blocks corresponding to qualified diagonals:
    ! L(#J,{ab})
    ! ---------------------------------------------------------

    do jVec=1,NumV
      do iAB=1,nQual(iSym)
        jAB = iQuAB(iAB,iSym)-iiBstR(iSym,2)
        W(jVec,iAB) = V(jAB,jVec+iVec0)
      end do
    end do

    ! Subtract:
    ! (gd|{ab}) <- (gd|{ab}) - sum_J L(gd,#J) * L(#J,{ab})
    ! for each ab in {ab}.
    ! ----------------------------------------------------

    call Cho_SubScr_Dia(V(:,iVec0+1),NumV,iSym,2,SSNorm)

    do iAB=1,nQual(iSym)
      do iShGD=1,nnShl
        nGD = nnBstRSh(iSym,iShGD,2)
        if (nGD < 1) cycle
        iGD = iiBstRSh(iSym,iShGD,2)
        xTot = xTot+1.0d0
        jAB = iQuab(iAB,iSym)-iiBstR(iSym,2)
        Tst = sqrt(DSPNm(iShGD)*DSubScr(jAB))
        if (Tst <= SSTau) cycle
        xDon = xDon+1.0d0
        call dGeMV_('N',nGD,NumV,xMOne,V(1+iGD:,iVec0+1),nnBstR(iSym,2),W(:,iAB),1,One,U(1+iGD:,iAB),1)
      end do
    end do

  else ! unscreened subtraction

    if (associated(LQ(iSym)%Array)) then

      ! If the qualified block, L({ab},#J), is already in core,
      ! use this block.
      ! -------------------------------------------------------

      call DGEMM_('N','T',nnBstR(iSym,2),nQual(iSym),NumV,xMOne,V(:,iVec0+1),nnBstR(iSym,2),LQ(iSym)%Array(:,iVec0+1), &
                  size(LQ(iSym)%Array,1),One,U,nnBstR(iSym,2))

    else

      lRow = nQual(iSym)
      lCol = NumV
      iS = 1
      iE = iS-1+lRow*lCol

      W(1:lRow,1:lCol) => WrK(iS:iE)

      ! Copy out sub-blocks corresponding to qualified diagonals:
      ! L({ab},#J).
      ! ---------------------------------------------------------

      do jVec=1,NumV
        do iAB=1,nQual(iSym)
          jAB = iQuAB(iAB,iSym)-iiBstR(iSym,2)
          W(iAB,jVec) = V(jAB,jVec+iVec0)
        end do
      end do

      ! Subtract:
      ! (gd|{ab}) <- (gd|{ab}) - sum_J L(gd,#J) * L({ab},#J)
      ! ----------------------------------------------------

      call DGEMM_('N','T',nnBstR(iSym,2),nQual(iSym),NumV,xMOne,V(:,iVec0+1),nnBstR(iSym,2),W,nQual(iSym),One,U,nnBstR(iSym,2))

    end if

  end if

  V => null()
  U => null()
  W => null()

end do

! Update statistics info.
! -----------------------

if (DoStat) nDGM_Call = nDGM_Call+nBatch
if (Cho_SScreen) then
  SubScrStat(1) = SubScrStat(1)+xTot
  SubScrStat(2) = SubScrStat(2)+xDon
end if

! Update global timing.
! ---------------------

if (DoTime) then
  call Cho_Timer(C2,W2)
  tDecom(1,3) = tDecom(1,3)+C2-C1
  tDecom(2,3) = tDecom(2,3)+W2-W1
end if

end subroutine Cho_VecBuf_Subtr
