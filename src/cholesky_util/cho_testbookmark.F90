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

use stdalloc

implicit none
integer irc
logical verbose
logical is1CCD
#include "cholesky.fh"
character*16 SecNam
parameter(SecNam='Cho_TestBookmark')
logical dealloc
integer nVec(8)
integer jrc
integer iSym
integer NumChoBak(8)
integer NumChTBak
real*8 delta(8)
real*8 Thr
real*8 ThrComBak
real*8 ErrMx
real*8, allocatable :: BkmDia(:), BkmDiaX(:)
logical Cho_1Center_Bak

irc = 0

if (verbose) call Cho_Head('Output from '//SecNam,'=',80,6)

! test 1: asking for accuracy below decomposition threshold
!         should result in failure (return code 1)
Thr = ThrCom*1.0d-1
call Cho_X_Bookmark(Thr,nSym,nVec,delta,jrc)
if (jrc == -1) then ! bookmarks not available
  irc = -1
  if (verbose) write(6,'(A)') 'Cho_X_Bookmark returned -1 [not available]. No further testing performed!'
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
Thr = -1.0d-12
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
if (ThrCom < 1.0d-4) then
  Thr = max(ThrCom*1.0d3,1.0d-14)
else
  Thr = max(ThrCom*1.0d2,1.0d-14)
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
  do iSym=1,nSym
    NumChoBak(iSym) = NumCho(iSym)
  end do
  NumChTBak = NumChT
  ThrComBak = ThrCom
  NumChT = 0
  ThrCom = 0.0d0
  do iSym=1,nSym
    NumCho(iSym) = nVec(iSym)
    NumChT = NumChT+NumCho(iSym)
    ThrCom = max(ThrCom,delta(iSym))
  end do
  call mma_allocate(BkmDia,nnBstRT(1),Label='BkmDia')
  call Cho_X_CalcChoDiag(jrc,BkmDia)
  if (jrc == 0) then
    call mma_allocate(BkmDiaX,nnBstRT(1),Label='BkmDiaX')
    call Cho_IODiag(BkmDiaX,2)
    call dAXPY_(nnBstRT(1),-1.0d0,BkmDia,1,BkmDiaX,1)
    if (Cho_1Center) then
      call Cho_TestBookmark_1CInit(dealloc)
      call Cho_MaxAbsDiag(BkmDiaX,1,ErrMx)
      if (dealloc) call Cho_TestBookmark_1CFinal()
    else
      call Cho_MaxAbsDiag(BkmDiaX,1,ErrMx)
    end if
    call FZero(DiaMax,nSym)
    call FZero(DiaMaxT,nSym)
    call mma_deallocate(BkmDiaX)
    if (abs(ErrMx-ThrCom) > 1.0d-12) then
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
  do iSym=1,nSym
    NumCho(iSym) = NumChoBak(iSym)
  end do
  NumChT = NumChTBak
  ThrCom = ThrComBak
  Cho_1Center = Cho_1Center_Bak
end if

end subroutine Cho_TestBookmark
