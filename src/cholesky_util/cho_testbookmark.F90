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

subroutine Cho_TestBookmark(irc,verbose,is1CCD)

use Cholesky, only: Cho_1Center, DiaMax, DiaMaxT, iAtomShl, nnBstRT, nSym, NumCho, NumChT, ThrCom
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
logical(kind=iwp), intent(in) :: verbose, is1CCD
logical(kind=iwp) :: Cho_1Center_Bak, dealloc
integer(kind=iwp) :: iSym, jrc, NumChoBak(8), NumChTBak, nVec(8)
real(kind=wp) :: delta(8), ErrMx, Thr, ThrComBak
real(kind=wp), allocatable :: BkmDia(:), BkmDiaX(:)
character(len=*), parameter :: SecNam = 'Cho_TestBookmark'

irc = 0

if (verbose) call Cho_Head('Output from '//SecNam,'=',80,u6)

! test 1: asking for accuracy below decomposition threshold
!         should result in failure (return code 1)
Thr = ThrCom*1.0e-1_wp
call Cho_X_Bookmark(Thr,nSym,nVec,delta,jrc)
if (jrc == -1) then ! bookmarks not available
  irc = -1
  if (verbose) write(u6,'(A)') 'Cho_X_Bookmark returned -1 [not available]. No further testing performed!'
  return
end if
if (jrc == 1) then
  if (verbose) call Cho_TestBookmark_Prt(1,'passed')
else
  irc = irc+1
  if (verbose) call Cho_TestBookmark_Prt(1,'failed')
end if

! test 2: asking for negative accuracy
!         should result in failure (return code 1)
Thr = -1.0e-12_wp
call Cho_X_Bookmark(Thr,nSym,nVec,delta,jrc)
if (jrc == 1) then
  if (verbose) call Cho_TestBookmark_Prt(2,'passed')
else
  irc = irc+1
  if (verbose) call Cho_TestBookmark_Prt(2,'failed')
end if

! test 3: asking for decomposition threshold should
!         give the total number of vectors.
Thr = ThrCom
call Cho_X_Bookmark(Thr,nSym,nVec,delta,jrc)
if (jrc == 0) then
  do iSym=1,nSym
    if (nVec(iSym) /= NumCho(iSym)) jrc = jrc+1
    if (delta(iSym) > ThrCom) jrc = jrc+1
  end do
  if (jrc == 0) then
    if (verbose) call Cho_TestBookmark_Prt(3,'passed')
  else
    irc = irc+1
    if (verbose) call Cho_TestBookmark_Prt(3,'failed')
  end if
else
  irc = irc+1
  if (verbose) call Cho_TestBookmark_Prt(3,'failed')
end if

! test 4: a threshold above the decomposition threshold
!         should result in fewer vectors and less accuracy.
if (ThrCom < 1.0e-4_wp) then
  Thr = max(ThrCom*1.0e3_wp,1.0e-14_wp)
else
  Thr = max(ThrCom*1.0e2_wp,1.0e-14_wp)
end if
call Cho_X_Bookmark(Thr,nSym,nVec,delta,jrc)
if (jrc == 0) then
  do iSym=1,nSym
    if (delta(iSym) > Thr) jrc = jrc+1
    if (nVec(iSym) > NumCho(iSym)) then
      jrc = jrc+1
    else if (nVec(iSym) == NumCho(iSym)) then
      if (delta(iSym) > ThrCom) jrc = jrc+1
    end if
  end do
  if (jrc == 0) then
    if (verbose) call Cho_TestBookmark_Prt(4,'passed')
  else
    irc = irc+1
    if (verbose) call Cho_TestBookmark_Prt(4,'failed')
  end if
else
  irc = irc+1
  if (verbose) call Cho_TestBookmark_Prt(4,'failed')
end if

! Test 5: check diagonal (only if previous tests passed).
if (irc /= 0) then
  if (verbose) call Cho_TestBookmark_Prt(5,'not executed')
else
  Cho_1Center_Bak = Cho_1Center
  Cho_1Center = is1CCD
  NumChoBak(1:nSym) = NumCho(1:nSym)
  NumChTBak = NumChT
  ThrComBak = ThrCom
  NumCho(1:nSym) = nVec(1:nSym)
  NumChT = sum(NumCho(1:nSym))
  ThrCom = Zero
  do iSym=1,nSym
    ThrCom = max(ThrCom,delta(iSym))
  end do
  call mma_allocate(BkmDia,nnBstRT(1),Label='BkmDia')
  call Cho_X_CalcChoDiag(jrc,BkmDia)
  if (jrc == 0) then
    call mma_allocate(BkmDiaX,nnBstRT(1),Label='BkmDiaX')
    call Cho_IODiag(BkmDiaX,2)
    BkmDiaX(:) = BkmDiaX(:)-BkmDia(:)
    if (Cho_1Center) then
      call Cho_TestBookmark_1CInit(dealloc)
      call Cho_MaxAbsDiag(BkmDiaX,1,ErrMx)
      if (dealloc .and. allocated(iAtomShl)) call mma_deallocate(iAtomShl)
    else
      call Cho_MaxAbsDiag(BkmDiaX,1,ErrMx)
    end if
    DiaMax(1:nSym) = Zero
    DiaMaxT(1:nSym) = Zero
    call mma_deallocate(BkmDiaX)
    if (abs(ErrMx-ThrCom) > 1.0e-12_wp) then
      irc = irc+1
      if (verbose) call Cho_TestBookmark_Prt(5,'failed')
    else
      if (verbose) call Cho_TestBookmark_Prt(5,'passed')
    end if
  else
    irc = irc+1
    if (verbose) call Cho_TestBookmark_Prt(5,'failed')
  end if
  call mma_deallocate(BkmDia)
  NumCho(1:nSym) = NumChoBak(1:nSym)
  NumChT = NumChTBak
  ThrCom = ThrComBak
  Cho_1Center = Cho_1Center_Bak
end if

end subroutine Cho_TestBookmark
