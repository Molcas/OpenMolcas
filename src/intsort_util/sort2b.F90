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
!***********************************************************************

subroutine SORT2B(iBin,nInts,iOrd,lSrtA,SrtArr,IOStk,lStk,nStk)
!***********************************************************************
!                                                                      *
!     Purpose: Store the sorted integrals on disk.                     *
!              Every record that is written contains a header of       *
!              length lTop which contains the number of integrals      *
!              on the record and its ordering number.                  *
!                                                                      *
!     Called from: Sort2                                               *
!                                                                      *
!     Calls to : PkI4,PkR8,DaFile                                      *
!                                                                      *
!     Calling parameters:                                              *
!     iBin   : Slice number                                            *
!     nInts  : no. of 2el integrals in slice                           *
!     iOrd   : ordering number of record                               *
!     SrtArr : Work space to keep the 2el integrals                    *
!                                                                      *
!     local data declarations:                                         *
!     PkVal  : I/O buffer contains packed integral values              *
!     IntLen : this buffer contains the length of the integrals        *
!              to be pack into the next buffer                         *
!                                                                      *
!*** M. Fuelscher and P.-Aa. Malmqvist, Univ. of Lund, Sweden, 1991 ****

use TwoDat, only: lStRec, lTop
use sort_data, only: iDaTwo, iDVBin, LuTwo, mDaTwo, nRec
use Constants, only: Zero
use Definitions, only: wp, iwp, u6, RtoB

implicit none
integer(kind=iwp), intent(in) :: iBin, nInts, lSrtA, lStk
integer(kind=iwp), intent(inout) :: iOrd, IOStk(lStk), nStk
real(kind=wp), intent(in) :: SrtArr(lSrtA)
integer(kind=iwp) :: iEnd, IntLen(4*lStRec), iOpt, iOptIO, iSave, iStart, iStk, jStk, kStk, llVRec, lVRec, mxVRec, nSave, nSaved
real(kind=wp) :: PkVal(lStRec), Dum(1)
#ifdef _DEBUGPRINT_
#include "print.fh"
integer(kind=iwp) :: iPrint, iRout
#endif

!----------------------------------------------------------------------*
!     pick up the print level                                          *
!----------------------------------------------------------------------*

#ifdef _DEBUGPRINT_
iRout = 86
iPrint = nPrint(iRout)
if (iPrint > 5) then
  write(u6,*) ' >>> Enter SORT2B <<<'
  write(u6,*) ' iBin  ',iBin
  write(u6,*) ' lSrtA ',lSrtA
  write(u6,*) ' nInts ',nInts
end if
if (iPrint >= 10) then
  call iVcPrt('stack of free records',' ',IOStk,nStk)
end if
#endif

!----------------------------------------------------------------------*
!     Transfer integrals from the sorting area to I/O buffer           *
!----------------------------------------------------------------------*
nRec(iBin) = 0
kStk = 0
nSaved = 0
iOpt = iDVBin(4,iBin)
do while (nSaved < nInts)
  iStart = nSaved+1
  iEnd = min(nSaved+4*lStRec,nInts)
  call R8Len(iOpt,1+iEnd-iStart,SrtArr(iStart),IntLen)
  mxVRec = (lStRec-lTop-1)*RtoB
  nSave = 0
  lVRec = 0
  do iSave=1,1+iEnd-iStart
    lVRec = lVRec+IntLen(iSave)
    if (lVRec < mxVRec) nSave = iSave
  end do
  if (nSave == 0) then
    write(u6,*)
    write(u6,'(2X,A,I3.3,A)') '*** Error in SORT2B ***'
    write(u6,'(2X,A)') 'nSave = 0'
    write(u6,*)
    call xFlush(u6)
    call Abend()
  end if
  lVRec = 0
  do iSave=1,nSave
    lVRec = lVRec+IntLen(iSave)
  end do
  if (lVRec > mxVRec) then
    write(u6,*)
    write(u6,'(2X,A,I3.3,A)') '*** Error in SORT2B ***'
    write(u6,'(2X,A)') 'An inconsistency has been deteced'
    write(u6,'(2X,A)') 'lVRec > mxVRec '
    write(u6,*)
    call xFlush(u6)
    call Abend()
  end if
  !--------------------------------------------------------------------*
  !   Pack integrals                                                   *
  !--------------------------------------------------------------------*
  call PKR8(iOpt,nSave,llVRec,SrtArr(iStart),PkVal(lTop+1))
  if (llVRec /= lVRec) then
    write(u6,*)
    write(u6,'(2X,A,I3.3,A)') '*** Error in SORT2B ***'
    write(u6,'(2X,A)') 'An inconsistency has been deteced'
    write(u6,'(2X,A)') 'llVBin # lVRec'
    write(u6,*)
    call xFlush(u6)
    call Abend()
  end if
  iOrd = iOrd+1
  PkVal(1) = 0
  PkVal(2) = iOrd
  PkVal(3) = nSave
  PkVal(4) = iOpt
  !--------------------------------------------------------------------*
  !   Write out the new record                                         *
  !--------------------------------------------------------------------*
  if (kStk < nStk) then
    ! use an old record, get address from IOStk
    kStk = kStk+1
    iDaTwo = IOStk(kStk)
  else
    iDaTwo = mDaTwo
    iOptIO = 0
    Dum(1) = Zero
    call dDAFILE(LuTwo,iOptIO,Dum,lStRec,mDaTwo)
  end if
# ifdef _DEBUGPRINT_
  if (iPrint >= 10) then
    write(u6,*) ' write record: iOrd,iDaTwo ',iOrd,iDaTwo
  end if
# endif

  ! Write out the buffer

  iOptIO = 1
  call dDAFILE(LuTwo,iOptIO,PkVal,lStRec,iDaTwo)

  nRec(iBin) = nRec(iBin)+1
  nSaved = nSaved+nSave

end do
!----------------------------------------------------------------------*
!     cleanup the stack of Record addresses                            *
!----------------------------------------------------------------------*
jStk = 0

! Do this only if the number of records was reduced

if (kStk < nStk) then

  ! Move unused record addresses to the start of IOStk!

  do iStk=kStk+1,nStk
    jStk = jStk+1
    IOStk(jStk) = IOStk(iStk)
  end do
end if
nStk = jStk

!----------------------------------------------------------------------*
!     exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine SORT2B
