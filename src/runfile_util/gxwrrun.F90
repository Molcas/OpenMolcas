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
! Copyright (C) 2003, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine writes a record into the runfile.                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************

subroutine gxWrRun(iRc,Label,data,nData,iOpt,RecTyp)

#include "runinfo.fh"
#include "runtypes.fh"
!----------------------------------------------------------------------*
! Declare arguments                                                    *
!----------------------------------------------------------------------*
integer iRc
character*(*) Label
character data(*)
integer nData
integer iOpt
integer RecTyp
!----------------------------------------------------------------------*
! Declare local data                                                   *
!----------------------------------------------------------------------*
character*64 ErrMsg
integer Lu
integer iDisk
integer DataAdr
logical ok
logical remove
integer item
integer i
integer NewLen

!----------------------------------------------------------------------*
! Check that arguments are ok.                                         *
!----------------------------------------------------------------------*
DataAdr = -99999
ok = .false.
if (RecTyp == TypInt) ok = .true.
if (RecTyp == TypDbl) ok = .true.
if (RecTyp == TypStr) ok = .true.
if (RecTyp == TypLgl) ok = .true.
if (.not. ok) call SysAbendMsg('gxWrRun','Argument RecTyp is of wrong type','Aborting')
if (nData < 0) call SysAbendMsg('gxWrRun','Number of data items less than zero','Aborting')
if (iOpt /= 0) then
  write(ErrMsg,*) 'Illegal option flag:',iOpt
  call SysAbendMsg('gxWrRun',ErrMsg,' ')
end if
iRc = 0
!----------------------------------------------------------------------*
! Does the runfile exist? If not create it.                            *
!----------------------------------------------------------------------*
call f_inquire(RunName,ok)
if (.not. ok) call MkRun(iRc,iOpt)
!----------------------------------------------------------------------*
! Open runfile.                                                        *
!----------------------------------------------------------------------*
call OpnRun(iRc,Lu,iOpt)
!----------------------------------------------------------------------*
! Do we have space left on file?                                       *
!----------------------------------------------------------------------*
if (RunHdr(ipItems) >= nToc) then
  call DaClos(Lu)
  call SysFilemsg('gxWrRun','Ran out of ToC record in RunFile',Lu,' ')
  call Abend()
end if
!----------------------------------------------------------------------*
! Read the ToC                                                         *
!----------------------------------------------------------------------*
iDisk = RunHdr(ipDaLab)
call cDaFile(Lu,icRd,TocLab,16*nToc,iDisk)
iDisk = RunHdr(ipDaPtr)
call iDaFile(Lu,icRd,TocPtr,nToc,iDisk)
iDisk = RunHdr(ipDaLen)
call iDaFile(Lu,icRd,TocLen,nToc,iDisk)
iDisk = RunHdr(ipDaMaxLen)
call iDaFile(Lu,icRd,TocMaxLen,nToc,iDisk)
iDisk = RunHdr(ipDaTyp)
call iDaFile(Lu,icRd,TocTyp,nToc,iDisk)
!----------------------------------------------------------------------*
! Reuse old field?                                                     *
!----------------------------------------------------------------------*
item = -1
do i=1,nToc
  if (TocLab(i) == Label) item = i
end do
NewLen = 0
if (item /= -1) then
  remove = .false.
  if (TocTyp(item) /= RecTyp) remove = .true.
  if (TocMaxLen(item) < nData) remove = .true.
  if (remove) then
    !write (6,*) '*******************************************'
    !write (6,'(a,a,a,i10)') 'Label=',Label,' expands in RUNFILE with size=',nData
    !write (6,*) '*******************************************'
    !call Abend()
    TocLab(item) = 'Empty   '
    TocPtr(item) = NulPtr
    TocLen(item) = 0
    TocTyp(item) = TypUnk
    item = -1
  else
    DataAdr = TocPtr(item)
    NewLen = TocLen(item)
  end if
  RunHdr(ipItems) = RunHdr(ipItems)-1
end if
!----------------------------------------------------------------------*
! Use new field?                                                       *
!----------------------------------------------------------------------*
if (item == -1) then
  do i=nToc,1,-1
    if (TocPtr(i) == NulPtr) item = i
  end do
  if (item == -1) then
    call DaClos(Lu)
    call SysFilemsg('gxWrRun','Internal inconsistency handling RunFile',Lu,' ')
    call Abend()
  end if
  DataAdr = RunHdr(ipNext)
end if
!----------------------------------------------------------------------*
! Write data to runfile and update header.                             *
!----------------------------------------------------------------------*
RunHdr(ipItems) = RunHdr(ipItems)+1
TocLab(item) = Label
TocPtr(item) = DataAdr
TocLen(item) = nData
TocMaxLen(item) = max(NewLen,nData)
TocTyp(item) = RecTyp

iDisk = DataAdr
call gzRWRun(Lu,icWr,data,nData,iDisk,RecTyp)

if (iDisk > RunHdr(ipNext)) RunHdr(ipNext) = iDisk
iDisk = 0
call iDaFile(Lu,icWr,RunHdr,nHdrSz,iDisk)
iDisk = RunHdr(ipDaLab)
call cDaFile(Lu,icWr,TocLab,16*nToc,iDisk)
iDisk = RunHdr(ipDaPtr)
call iDaFile(Lu,icWr,TocPtr,nToc,iDisk)
iDisk = RunHdr(ipDaLen)
call iDaFile(Lu,icWr,TocLen,nToc,iDisk)
iDisk = RunHdr(ipDaMaxLen)
call iDaFile(Lu,icWr,TocMaxLen,nToc,iDisk)
iDisk = RunHdr(ipDaTyp)
call iDaFile(Lu,icWr,TocTyp,nToc,iDisk)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
call DaClos(Lu)

return

end subroutine gxWrRun
