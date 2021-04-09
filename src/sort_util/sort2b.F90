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
!              length lTop which containd the number of integrals      *
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
!     Global data declarations (Include files) :                       *
!     TwoDef  : definitions of the record structure                    *
!     Srt0    : common block containing information pertinent to       *
!               the calculation of 2el integral sequence numbers       *
!     Srt1    : common block containing information the number of      *
!               bins and partitioning of symmetry blocks               *
!     Srt2    : common block containing information pertinent to       *
!               the bin sorting algorithm                              *
!     Srt3    : dynamic stack to control inout and output of           *
!               integral buffers                                       *
!                                                                      *
!     local data declarations:                                         *
!     PkVal  : I/O buffer contains packed integral values              *
!     IntLen : this buffer contains the length of the integrals        *
!              to be pack into the next buffer                         *
!                                                                      *
!*** M. Fuelscher and P.-Aa. Malmqvist, Univ. of Lund, Sweden, 1991 ****

use srt2
implicit real*8(A-H,O-Z)

#include "TwoDat.fh"
#include "srt0.fh"
#include "srt1.fh"

#include "SysDef.fh"
#include "print.fh"

dimension PkVal(lStRec)
dimension IntLen(4*lStRec)
dimension SrtArr(lSrtA)
integer IOStk(lStk)

!----------------------------------------------------------------------*
!     pick up the print level                                          *
!----------------------------------------------------------------------*

#ifdef _DEBUGPRINT_
iRout = 86
iPrint = nPrint(iRout)
if (iPrint > 5) then
  write(6,*) ' >>> Enter SORT2B <<<'
  write(6,*) ' iBin  ',iBin
  write(6,*) ' LSrtA ',lSrtA
  write(6,*) ' nInts ',nInts
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
    write(6,*)
    write(6,'(2X,A,I3.3,A)') '*** Error in SORT2B ***'
    write(6,'(2X,A)') 'nSave = 0'
    write(6,*)
    call xFlush(6)
    call Abend()
  end if
  lVRec = 0
  do iSave=1,nSave
    lVRec = lVRec+IntLen(iSave)
  end do
  if (lVRec > mxVRec) then
    write(6,*)
    write(6,'(2X,A,I3.3,A)') '*** Error in SORT2B ***'
    write(6,'(2X,A)') 'An inconsistency has been deteced'
    write(6,'(2X,A)') 'lVRec > mxVRec '
    write(6,*)
    call xFlush(6)
    call Abend()
  end if
  !--------------------------------------------------------------------*
  !   Pack integrals                                                   *
  !--------------------------------------------------------------------*
  call PKR8(iOpt,nSave,llVRec,SrtArr(iStart),PkVal(lTop+1))
  if (llVRec /= lVRec) then
    write(6,*)
    write(6,'(2X,A,I3.3,A)') '*** Error in SORT2B ***'
    write(6,'(2X,A)') 'An inconsistency has been deteced'
    write(6,'(2X,A)') 'llVBin # lVRec'
    write(6,*)
    call xFlush(6)
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
    call dDAFILE(LuTwo,iOptIO,[0.0d0],lStRec,mDaTwo)
  end if
# ifdef _DEBUGPRINT_
  if (iPrint >= 10) then
    write(6,*) ' write record: iOrd,iDaTwo ',iOrd,iDaTwo
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
