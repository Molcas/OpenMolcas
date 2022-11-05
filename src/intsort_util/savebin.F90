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
! Copyright (C) 1991, Markus P. Fuelscher                              *
!               1991, Per Ake Malmqvist                                *
!               1998, Roland Lindh                                     *
!***********************************************************************

subroutine SaveBin(iBin,iOpt)
!***********************************************************************
!                                                                      *
!     Purpose: Phase 1 of the bin sorting algorithm                    *
!              Once a bin is filled it has to be dumped to disk.       *
!              First, however, we like to sort and pack the buffers.   *
!              A further complication arises due to the fact that a    *
!              bin is shorter than the records as used in the wave-    *
!              function codes.                                         *
!                                                                      *
!     Called from: SORT1A and SORT1B                                   *
!                                                                      *
!     Calls to : PKI4,PKR8,DaFile,SetVec,ISORTX,I4Len,R8Len            *
!                                                                      *
!     Calling Parameters:                                              *
!     iBin   : Bin number to be saved                                  *
!                                                                      *
!     local data declarations:                                         *
!     Scr    : temporary array used for initializing records           *
!                                                                      *
!     Modified to remove sorting step, R. Lindh, March '98             *
!                                                                      *
!*** M. Fuelscher and P.-Aa. Malmqvist, Univ. of Lund, Sweden, 1991 ****

use sort_data, only: iDaTmp, iDATwo, iDIBin, iDVBin, IndBin, lIndx, lInts, LuTmp, LuTwo, lwIBin, lwVBin, mDaTmp, mDaTwo, mInt, &
                     nRec, n_Int, ValBin
use TwoDat, only: lDaRec, lStRec, lTop, nSect
use Pack_mod, only: isPack
use Constants, only: Zero
use Definitions, only: wp, iwp, u6, ItoB, RtoB

implicit none
integer(kind=iwp), intent(in) :: iBin, iOpt
integer(kind=iwp) :: i, idiv, iOptIO, iSave, iScr(lStRec), lIBin, lIRec, lVBin, lVRec, mInds, mInts, mxIRec, mxVRec, nInts, nKeep, &
                     nSave
real(kind=wp) :: Scr(lStRec)

!----------------------------------------------------------------------*
!         as the packed integral labels add an extra 1-2 Byte          *
!         disk space per integral we have to adjust the record         *
!         length of LuTmp to the different machines.                   *
!----------------------------------------------------------------------*
idiv = ItoB/2
if (isPack) idiv = idiv/2
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
nInts = n_Int(iBin)
!----------------------------------------------------------------------*
!         precompute the packed length of a record                     *
!----------------------------------------------------------------------*
call I4LEN(nInts,lwIBin(1,iBin),lIndx)
call R8LEN(iOpt,nInts,lwVBin(1,iBin),lInts)
mxIRec = ((lDaRec/idiv)-lTop-1)*ItoB
mxVRec = (lDaRec-lTop-1)*RtoB
nSave = 0
lIRec = 0
lVRec = 0
do iSave=1,nInts
  lIRec = lIRec+lIndx(iSave)
  lVRec = lVRec+lInts(iSave)
  if ((lIRec < mxIRec) .and. (lVRec < mxVRec)) nSave = iSave
end do
lIRec = 0
lVRec = 0
do iSave=1,nSave
  lIRec = lIRec+lIndx(iSave)
  lVRec = lVRec+lInts(iSave)
end do
!----------------------------------------------------------------------*
!         now pack and check the packed record length again            *
!----------------------------------------------------------------------*
call PKI4(nSave,lIBin,lwIBin(1,iBin),IndBin(lTop+1))
mInds = (lIBin+ItoB-1)/ItoB
mInt(3,iBin) = mInt(3,iBin)+mInds
if (lIBin > mxIRec) then
  write(u6,*)
  write(u6,'(2X,A,I3.3,A)') '*** Error in SAVEBIN ***'
  write(u6,'(2X,A)') 'An inconsistency has been deteced'
  write(u6,'(2X,A)') 'lIRec > mxIRec '
  write(u6,*)
  call xFlush(u6)
  call Abend()
end if
if (lIBin /= lIRec) then
  write(u6,*)
  write(u6,'(2X,A,I3.3,A)') '*** Error in SAVEBIN ***'
  write(u6,'(2X,A)') 'An inconsistency has been deteced'
  write(u6,'(2X,A)') 'lIBin # lIRec'
  write(u6,*)
  call xFlush(u6)
  call Abend()
end if
call PKR8(iOpt,nSave,lVBin,lwVBin(1,iBin),ValBin(lTop+1))
mInts = (lVBin+RtoB-1)/RtoB
mInt(2,iBin) = mInt(2,iBin)+mInts
if (lVBin > mxVRec) then
  write(u6,*)
  write(u6,'(2X,A,I3.3,A)') '*** Error in SAVEBIN ***'
  write(u6,'(2X,A)') 'An inconsistency has been deteced'
  write(u6,'(2X,A)') 'lVRec > mxVRec '
  write(u6,*)
  call xFlush(u6)
  call Abend()
end if
if (lVBin /= lVRec) then
  write(u6,*)
  write(u6,'(2X,A,I3.3,A)') '*** Error in SAVEBIN ***'
  write(u6,'(2X,A)') 'An inconsistency has been deteced'
  write(u6,'(2X,A)') 'lVBin # lVRec'
  write(u6,*)
  call xFlush(u6)
  call Abend()
end if
!----------------------------------------------------------------------*
!     write the record header                                          *
!----------------------------------------------------------------------*
IndBin(1) = iDIBin(2,iBin)
IndBin(2) = lIBin+lTop*ItoB
IndBin(3) = nSave
ValBin(1) = iDVBin(2,iBin)
ValBin(2) = lVBin+lTop*RtoB
ValBin(3) = nSave
!----------------------------------------------------------------------*
!     save the record on disk                                          *
!----------------------------------------------------------------------*
iDaTmp = iDIBin(1,iBin)
iDaTwo = iDVBin(1,iBin)
iOptIO = 1

! If a new sector mark it on the disk to make sure that it is continuous.

if (mod(nRec(iBin),nSect) == 0) then
  iDaTmp = mDaTmp
  iDaTwo = mDaTwo

  iScr(:) = 0
  call iDAFILE(LuTmp,iOptIO,iScr,(lStRec/idiv),mDaTmp)

  Scr(:) = Zero
  call dDAFILE(LuTwo,iOptIO,Scr,lStRec,mDaTwo)

  iDIBin(2,iBin) = iDaTmp
  iDVBin(2,iBin) = iDaTwo
  if (iDVBin(3,iBin) < 0) iDVBin(3,iBin) = iDaTwo
end if

call iDAFILE(LuTmp,iOptIO,IndBin,(lDaRec/idiv),iDaTmp)
iDIBin(1,iBin) = iDaTmp
call dDAFILE(LuTwo,iOptIO,ValBin,lDaRec,iDaTwo)
iDVBin(1,iBin) = iDaTwo
nRec(iBin) = nRec(iBin)+1
!----------------------------------------------------------------------*
!         cleanup                                                      *
!----------------------------------------------------------------------*
nKeep = nInts-nSave
if (nKeep > 0) then
  do i=1,nKeep
    lwIBin(i,iBin) = lwIBin(nSave+i,iBin)
    lwVBin(i,iBin) = lwVBin(nSave+i,iBin)
  end do
end if
n_Int(iBin) = nKeep

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine SaveBin
