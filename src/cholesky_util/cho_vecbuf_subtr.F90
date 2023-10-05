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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine Cho_VecBuf_Subtr(xInt,Wrk,lWrk,iSym,DoTime,DoStat)
!
! Purpose: subtract contributions to qualified columns from the
!          vectors stored in the buffer (if any).
!
! DoTime: time as vector subtraction.
! DpStat: update statistics info (#calls to dGeMM).

use Cholesky, only: Cho_SScreen, CHVBUF, DSPNm, DSubScr, iiBstR, iiBstRSh, ip_CHVBUF_SYM, iQuAB, l_CHVBUF_SYM, LQ, LuPri, &
                    nDGM_call, nnBstR, nnBstRSh, nnShl, nQual, nVec_in_Buf, SSNorm, SSTau, SubScrStat, TDECOM
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), target, intent(inout) :: xInt(*)
integer(kind=iwp), intent(in) :: lWrk, iSym
real(kind=wp), target, intent(out) :: Wrk(lWrk)
logical(kind=iwp), intent(in) :: DoTime, DoStat
integer(kind=iwp) :: iAB, iBatch, iE, iGD, iS, iShGD, iVec0, jAB, jVec, lCol, lRow, nBatch, nGD, NumV, nVec
real(kind=wp) :: C1, C2, Tst, W1, W2, xDon, xTot
real(kind=wp), pointer :: U(:,:), V(:,:), W(:,:)
#ifdef _DEBUGPRINT_
#define _DBG_ .true.
#else
#define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: LocDbg = _DBG_
character(len=*), parameter :: SecNam = 'Cho_VecBuf_Subtr'

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

if (DoTime) call CWTime(C1,W1)

! Initialize.
! -----------

xTot = Zero
xDon = Zero

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
        xTot = xTot+One
        jAB = iQuab(iAB,iSym)-iiBstR(iSym,2)
        Tst = sqrt(DSPNm(iShGD)*DSubScr(jAB))
        if (Tst <= SSTau) cycle
        xDon = xDon+One
        call dGeMV_('N',nGD,NumV,-One,V(1+iGD:,iVec0+1),nnBstR(iSym,2),W(:,iAB),1,One,U(1+iGD:,iAB),1)
      end do
    end do

  else ! unscreened subtraction

    if (associated(LQ(iSym)%A)) then

      ! If the qualified block, L({ab},#J), is already in core,
      ! use this block.
      ! -------------------------------------------------------

      call DGEMM_('N','T',nnBstR(iSym,2),nQual(iSym),NumV,-One,V(:,iVec0+1),nnBstR(iSym,2),LQ(iSym)%A(:,iVec0+1), &
                  size(LQ(iSym)%A,1),One,U,nnBstR(iSym,2))

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

      call DGEMM_('N','T',nnBstR(iSym,2),nQual(iSym),NumV,-One,V(:,iVec0+1),nnBstR(iSym,2),W,nQual(iSym),One,U,nnBstR(iSym,2))

    end if

  end if

  nullify(V)
  nullify(W)

end do

nullify(U)

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
  call CWTime(C2,W2)
  tDecom(1,3) = tDecom(1,3)+C2-C1
  tDecom(2,3) = tDecom(2,3)+W2-W1
end if

end subroutine Cho_VecBuf_Subtr
