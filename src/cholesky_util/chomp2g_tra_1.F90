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

subroutine ChoMP2g_Tra_1(COrb1,COrb2,Diag,DoDiag,Wrk,lWrk,iSym,iMoType1,iMoType2)
!
! Thomas Bondo Pedersen, Dec. 2010.
!
! Purpose: transform Cholesky vectors to (pq) MO basis for symmetry
!          block iSym. Files are assumed open.
!          If requested (DoDiag=.true.), compute (pq|pq) integral
!          diagonal.

use Cholesky, only: InfVec, nnBstR, NumCho
use ChoMP2, only: iAdrOff, lUnit_F, nAdrOff, nMoAo, nMoMo, nMoType
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: COrb1(*), COrb2(*)
integer(kind=iwp), intent(in) :: lWrk, iSym, iMoType1, iMoType2
real(kind=wp), intent(_OUT_) :: Diag(*)
logical(kind=iwp), intent(in) :: DoDiag
real(kind=wp), intent(out) :: Wrk(lWrk)
integer(kind=iwp) :: iAdr, iBat, iLoc, iOpt, irc, iRed, iRedC, iVec, iVec1, iVec2, iVecType, jNum, jVec, jVec1, kChoAO, kChoMO, &
                     kEnd0, kHlfTr, kOff, kOffMO, lChoAO, lChoMO, lHlfTr, lRead, lWrk0, lWrk1, mUsed, nMOVec, NumBat, NumV
character(len=*), parameter :: SecNam = 'ChoMP2_Tra_1'
integer(kind=iwp), external :: Cho_lRead

! Check what type of Cholesky vector to make (fro-occ, occ-occ.....)
iVecType = iMoType2+(iMoType1-1)*nMoType

! Check if anything to do.
! ------------------------

if ((NumCho(iSym) < 1) .or. (nMoMo(iSym,iVecType) < 1)) return

! Initialize Diag (if needed).
! ----------------------------

if (DoDiag) Diag(1:nMoMo(iSym,iVecType)) = Zero

! Allocate memory for half-transformed vector.
! --------------------------------------------

lHlfTr = nMoAo(iSym,iMoType1)

kHlfTr = 1
kEnd0 = kHlfTr+lHlfTr
lWrk0 = lWrk-kEnd0+1
if (lWrk0 < (nMoMo(iSym,iVecType)+nnBstR(iSym,1))) call SysAbendMsg(SecNam,'insufficient memory','[0]')

! Reserve memory for reading AO vectors.
! --------------------------------------

lRead = Cho_lRead(iSym,lWrk0)
if (lRead < 1) then
  write(u6,*) SecNam,': memory error: lRead = ',lRead
  call SysAbendMsg(SecNam,'memory error',' ')
  lWrk1 = 0 ! to avoid compiler warnings...
else
  lWrk1 = lWrk0-lRead
  if (lWrk1 < nMoMo(iSym,iVecType)) then
    lWrk1 = nMoMo(iSym,iVecType)
    lRead = lWrk-nMoMo(iSym,iVecType)
  end if
end if

! Set up batch.
! -------------

nMOVec = min(lWrk1/nMoMo(iSym,iVecType),NumCho(iSym))
if (nMOVec < 1) call SysAbendMsg(SecNam,'insufficient memory','[1]')
NumBat = (NumCho(iSym)-1)/nMOVec+1

! Set reduced set handles.
! ------------------------

iRedC = -1
iLoc = 3

! Transform each batch of vectors and compute diagonal contributions
! (if requested).
! ------------------------------------------------------------------

do iBat=1,NumBat

  if (iBat == NumBat) then
    NumV = NumCho(iSym)-nMOVec*(NumBat-1)
  else
    NumV = nMOVec
  end if
  iVec1 = nMOVec*(iBat-1)+1
  iVec2 = iVec1+NumV-1

  lChoMO = nMoMo(iSym,iVecType)*NumV

  kChoMO = kEnd0
  kChoAO = kChoMO+lChoMO
  lChoAO = lWrk0-kChoAO+1

  kOffMO = kChoMO
  jVec1 = iVec1
  do while (jVec1 <= iVec2)

    jNum = 0
    call Cho_VecRd(Wrk(kChoAO),lChoAO,jVec1,iVec2,iSym,jNum,iRedC,mUsed)
    if (jNum < 1) call SysAbendMsg(SecNam,'insufficient memory','[2]')

    kOff = kChoAO
    do jVec=1,jNum
      iVec = jVec1+jVec-1
      iRed = InfVec(iVec,2,iSym)
      if (iRedC /= iRed) then
        irc = 0
        call Cho_X_SetRed(irc,iLoc,iRed)
        if (irc /= 0) call SysAbendMsg(SecNam,'error in Cho_X_SetRed',' ')
        iRedC = iRed
      end if

      call ChoMP2g_TraVec(Wrk(kOff),Wrk(kOffMO),COrb1,COrb2,Wrk(kHlfTr),lHlfTr,iSym,1,1,iLoc,iMoType1,iMoType2)
      kOff = kOff+nnBstR(iSym,iLoc)
      kOffMO = kOffMO+nMoMo(iSym,iVecType)
    end do

    jVec1 = jVec1+jNum

  end do

  iOpt = 1
  iAdr = nAdrOff(iSym)+nMoMo(iSym,iVecType)*(iVec1-1)+1
  iAdrOff(iSym,iVecType) = nAdrOff(iSym)
  call ddaFile(lUnit_F(iSym,1),iOpt,Wrk(kChoMO),lChoMO,iAdr)

  if (DoDiag) then
    do iVec=1,NumV
      kOff = kChoMO+nMoMo(iSym,iVecType)*(iVec-1)-1
      Diag(1:nMoMo(iSym,iVecType)) = Diag(1:nMoMo(iSym,iVecType))+Wrk(kOff+1:kOff+nMoMo(iSym,iVecType))**2
    end do
  end if

end do
! When we reach this point we have written all vectors of the present
! type for this symmetry and need to remember were we should continue
! to write the next type.
if (iVecType /= 9) nAdrOff(iSym) = iAdr-1

end subroutine ChoMP2g_Tra_1
