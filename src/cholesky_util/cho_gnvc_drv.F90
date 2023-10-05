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

subroutine Cho_GnVc_Drv(irc,Diag)
!
! Purpose: generate vectors from a known "map" and diagonal.
!          First reduced set must have been set up, which is
!          reasonable, since it is naturally done along with the
!          diagonal.

use Data_Structures, only: Alloc1DiArray_Type
use Cholesky, only: DID_DECDRV, iiBstR, INF_PASS, InfVec, iOff_col, IPRINT, iQuAB, LuPri, nnBstR, nnShl, nQual, nSym, NumCho, &
                    NumChT, TDECDRV, XnPass
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(inout) :: Diag(*)
integer(kind=iwp) :: i, iAB, iPass, iPass1, iPass2, iSym, iV, iV1, iV2, iVec1, jAB, jPass, jRed, kAB, l_Int, l_Wrk, l_WrkT, &
                     LastRed(8), lThis, nBatch, nPass, nScrV(8), nTotVec, NumInt, NumPass, NumSP, nVec
real(kind=wp) :: dl_Int, dl_WrkT, tCPU1, tCPU2, TlDec, TlDec1, TlDec2, TlInt, TlInt1, TlInt2, TlTot, TlTot1, TlTot2, tWall1, &
                 tWall2, WlDec, WlDec1, WlDec2, WlInt, WlInt1, WlInt2, WlTot, WlTot1, WlTot2
character(len=26) :: String
character(len=2) :: Unt
type(Alloc1DiArray_Type) :: RS2RS(8)
integer(kind=iwp), allocatable :: iVecRS(:,:), LISTSP(:), nVecRS(:,:)
real(kind=wp), allocatable :: Wrk(:), xInt(:)
character(len=*), parameter :: SecNam = 'Cho_GnVc_Drv'

! Start timing.
! -------------

call CWTime(tCPU1,tWall1)

! Set return code.
! ----------------

irc = 0

! Copy first reduced set to current reduced set stored at location 2.
! -------------------------------------------------------------------

call Cho_X_RSCopy(irc,1,2)
if (irc /= 0) then
  write(Lupri,*) SecNam,': Cho_X_RSCopy returned ',irc
  call Finish_this()
  return
end if

! Allocate memory for shell pair list.
! ------------------------------------

call mma_allocate(ListSP,nnShl,Label='ListSP')

! Set number of integral passes (= number of reduced sets).
! Set ID of last reduced set in each symmetry.
! ---------------------------------------------------------

nPass = 0
do iSym=1,nSym
  if (NumCho(iSym) < 1) then
    LastRed(iSym) = 0
  else
    LastRed(iSym) = InfVec(NumCho(iSym),2,iSym)
  end if
  nPass = max(nPass,LastRed(iSym))
end do
if (nPass < 1) then
  call Cho_Quit('nPass is non-positive in '//SecNam,102)
else if (nPass /= XnPass) then
  call Cho_Quit('nPass != XnPass in '//SecNam,102)
end if

! nVecRS(iSym,jRed): #vectors of sym. iSym, reduced set jRed.
! iVecRS(iSym,jRed): 1st vec. of sym. iSym, reduced set jRed.
!                        (0 means no vectors!!)
! -----------------------------------------------------------

call mma_allocate(nVecRS,nSym,nPass,Label='nVecRS')
call mma_allocate(iVecRS,nSym,nPass,Label='iVecRS')
nVecRS(:,:) = 0
iVecRS(:,:) = 0

do iSym=1,nSym
  nTotVec = 0
  do jRed=1,LastRed(iSym)
    call Cho_X_nVecRS(jRed,iSym,iVec1,nVec)
    if ((nVec < 0) .or. (iVec1 < 0)) then
      call Cho_Quit(SecNam//'Cho_X_nVecRS returned negative nVec',103)
    else
      nVecRS(iSym,jRed) = nVec
      iVecRS(iSym,jRed) = iVec1
    end if
    nTotVec = nTotVec+nVecRS(iSym,jRed)
  end do
  if (nTotVec /= NumCho(iSym)) then
    call Cho_Quit('Initialization error in '//SecNam,102)
  end if
end do

! Allocate RS-to-RS mapping array.
! --------------------------------

do iSym=1,nSym
  call mma_allocate(RS2RS(iSym)%A,nnBstR(iSym,1),Label='RS2RS(iSym)%A')
end do

! Split remaining memory.
! -----------------------

call mma_MaxDBLE(l_WrkT)
call mma_allocate(Wrk,l_WrkT,Label='Wrk')
if (l_WrkT < 2) then
  write(Lupri,*) SecNam,': max. allocatable memory is ',l_WrkT
  write(Lupri,*) 'Please increase available memory.'
  irc = 101
  call Finish_this()
  return
else
  l_Wrk = l_WrkT/2
end if

#ifdef _DEBUGPRINT_
! Debug: force batching.
! ----------------------

lNdMx = 0
do iPass=1,nPass
  lNeed = sum(nnBstR(1:nSym,2)*nVecRS(1:nSym,iPass))
  lNdMx = max(lNdMx,lNeed)
end do
l_Wrk = min(l_Wrk,lNdMx)
#endif

! Reinitialize vector counters.
! -----------------------------

NumCho(1:nSym) = 0
NumChT = 0

! Start batch loop over integral passes.
! --------------------------------------

iPass = 0
nBatch = 0
do while (iPass < nPass)

  if (iPrint >= INF_PASS) call CWTime(TlTot1,WlTot1)

  ! Update batch counter.
  ! ---------------------

  nBatch = nBatch+1
  iPass1 = iPass+1
  if (nBatch > nPass) then
    write(Lupri,*) SecNam,': batch counter exceeds pass counter: ',nBatch,' > ',nPass
    irc = 103
    call Finish_this()
    return
  end if

  ! Set up this batch based on current reduced set.
  ! NumPass: the number of passes treated in this batch.
  ! --------------------------------------------------------

  jRed = iPass
  NumPass = 0
  l_Int = 0
  lThis = 0
  do while (jRed < nPass)
    jRed = jRed+1
    lThis = sum(nnBstR(1:nSym,2)*nVecRS(1:nSym,jRed))
    l_Int = l_Int+lThis
    if (l_Int > l_Wrk) then
      l_Int = l_Int-lThis ! reset memory need
      jRed = nPass ! break loop
    else
      NumPass = NumPass+1
    end if
  end do
  if (NumPass < 1) then
    write(Lupri,*) SecNam,': insufficient memory for batch ',nBatch
    write(Lupri,*) SecNam,': at least  ',lThis,' double precision words needed.'
    write(Lupri,*) SecNam,': available ',l_Wrk,' double precision words.'
    irc = 101
    call Finish_this()
    return
  end if

  ! Print.
  ! ------

  if (iPrint >= INF_PASS) then
    write(String,'(A19,I7)') 'Integral Pass Batch',nBatch
    call Cho_Head(String,'*',80,Lupri)
    write(Lupri,'(A,I5,A,I5,A,I5,A)') 'Integral passes treated:',iPass1,' to',iPass+NumPass,' of',nPass,' passes.'
    call Cho_Word2Byte(l_WrkT,8,dl_WrkT,Unt)
    write(Lupri,'(A,I10,A,F10.3,A,A)') 'Total memory available           : ',l_WrkT,' 8-byte words; ',dl_WrkT,' ',Unt
    call Cho_Word2Byte(l_Int,8,dl_Int,Unt)
    write(Lupri,'(A,I10,A,F10.3,A,A)') 'Memory used for integrals/vectors: ',l_Int,' 8-byte words; ',dl_Int,' ',Unt
    nScrV(1:nSym) = 0
    do i=iPass1,iPass+NumPass
      nScrV(1:nSym) = nScrV(1:nSym)+nVecRS(1:nSym,i)
    end do
    write(Lupri,'(A,8I8)') '#vec. gener.  : ',(nScrV(i),i=1,nSym)
    call XFlush(Lupri)
  end if

  ! Allocate memory for integral columns and initialize.
  ! ----------------------------------------------------

  call mma_allocate(xInt,l_Int,Label='xInt')
  xInt(:) = Zero

  ! Set up map from first reduced set to reduced set iPass1.
  ! --------------------------------------------------------

  irc = 0
  call Cho_X_RSCopy(irc,1,3)
  if (irc /= 0) call Cho_Quit(SecNam//': non-zero return code from Cho_X_RSCopy',104)
  do iSym=1,nSym
    call Cho_RS2RS(RS2RS(iSym)%A,size(RS2RS(iSym)%A),3,2,iPass1,iSym)
  end do

  ! Set up "qualified column" index arrays.
  ! nQual(iSym): #qualifieds in symmetry iSym
  ! iQuAB(iAB,iSym): addr of qualified iAB, sym. iSym in curr.
  !                  reduced set.
  ! ----------------------------------------------------------

  nQual(1:nSym) = 0
  iPass2 = iPass1+NumPass-1
  do jPass=iPass1,iPass2
    do iSym=1,nSym
      iV1 = iVecRS(iSym,jPass)
      iV2 = iV1+nVecRS(iSym,jPass)-1
      do iV=iV1,iV2
        iAB = InfVec(iV,1,iSym) ! addr in 1st red. set
        jAB = iAB-iiBstR(iSym,1)
        kAB = RS2RS(jAB)%A(iSym) ! addr in curr. red. set
#       ifdef _DEBUGPRINT_
        if ((kAB < 1) .or. (kAB > nnBstR(iSym,2))) then
          write(Lupri,*) SecNam,': illegal kAB = ',kAB
          write(Lupri,*) 'Vector, symmetry, pass: ',IV,iSym,jPass
          write(Lupri,*) 'Allowed range: 1 - ',nnBstR(iSym,2)
          call Cho_Quit('Index error in '//SecNam,104)
        end if
#       endif
        nQual(iSym) = nQual(iSym)+1
        iQuAB(nQual(iSym),iSym) = iiBstR(iSym,2)+kAB
      end do
    end do
  end do

  iOff_Col(1) = 0
  NumInt = nnBstR(1,2)*nQual(1)
  do iSym=2,nSym
    iOff_Col(iSym) = NumInt
    NumInt = NumInt+nnBstR(iSym,2)*nQual(iSym)
  end do
  if (l_Int /= NumInt) call Cho_Quit('Integral memory error in '//SecNam,101)

  ! Compute all integrals needed for NumPass passes.
  ! ------------------------------------------------

  if (iPrint >= INF_PASS) call CWTime(TlInt1,WlInt1)
  NumSP = 0
  call Cho_GnVc_GetInt(xInt,size(xInt),nVecRS,iVecRS,ListSP,nSym,nPass,nnShl,iPass1,NumPass,NumSP)
  if (NumSP < 1) call Cho_Quit('No shell pairs calculated!',104)
  if (iPrint >= INF_PASS) call CWTime(TlInt2,WlInt2)

  ! Generate vectors.
  ! -----------------

  if (iPrint >= INF_PASS) call CWTime(TlDec1,WlDec1)
  call Cho_GnVc_GenVec(Diag,xInt,size(xInt),nVecRS,iVecRS,RS2RS,nSym,nPass,iPass1,NumPass)
  if (iPrint >= INF_PASS) call CWTime(TlDec2,WlDec2)

  ! Deallocate memory.
  ! ------------------

  call mma_deallocate(xInt)

  ! Print timing for this batch.
  ! ----------------------------

  if (iPrint >= INF_PASS) then
    TlInt = TlInt2-TlInt1
    WlInt = WlInt2-WlInt1
    TlDec = TlDec2-TlDec1
    WlDec = WlDec2-WlDec1
    call CWTime(TlTot2,WlTot2)
    TlTot = TlTot2-TlTot1
    WlTot = WlTot2-WlTot1
    write(Lupri,'(/,A,I7,A)') 'Overall timings for integral pass batch',nBatch,' (CPU/Wall in seconds):'
    write(Lupri,'(A,F12.2,1X,F12.2)') 'Integrals    : ',TlInt,WlInt
    write(Lupri,'(A,F12.2,1X,F12.2)') 'Decomposition: ',TlDec,WlDec
    write(Lupri,'(A,F12.2,1X,F12.2)') 'Total        : ',TlTot,WlTot
    write(Lupri,'(A,I7,A,I7,A)') 'Integral passes treated:',iPass1,' to',iPass1-1+NumPass
    call XFlush(Lupri)
  end if

  ! Update pass counter.
  ! --------------------

  iPass = iPass+NumPass

end do ! integral pass loop

! Exit after deallocating memory.
! -------------------------------

call Finish_this()

return

contains

subroutine Finish_this()

  integer(kind=iwp) :: iSym

  Did_DecDrv = .true.

  call mma_deallocate(Wrk)
  do iSym=1,nSym
    call mma_deallocate(RS2RS(iSym)%A)
  end do
  call mma_deallocate(iVecRS)
  call mma_deallocate(nVecRS)
  call mma_deallocate(ListSP)

  call CWTime(tCPU2,tWall2)
  tDecDrv(1) = tDecDrv(1)+tCPU2-tCPU1
  tDecDrv(2) = tDecDrv(2)+tWall2-tWall1

end subroutine Finish_this

end subroutine Cho_GnVc_Drv
