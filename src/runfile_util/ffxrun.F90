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
! This routine locates a field in the runfile.                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************

subroutine ffxRun(iRc,Label,nData,RecTyp,iOpt)

#include "runinfo.fh"
#include "runtypes.fh"
#include "runrc.fh"
!----------------------------------------------------------------------*
! Declare arguments                                                    *
!----------------------------------------------------------------------*
integer iRc
character*(*) Label
integer nData
integer RecTyp
integer iOpt
!----------------------------------------------------------------------*
! Declare local data                                                   *
!----------------------------------------------------------------------*
character*64 ErrMsg
character*16 CmpLab1
character*16 CmpLab2
integer Lu
integer iDisk
logical ok
integer item
integer i

!----------------------------------------------------------------------*
! Check that arguments are ok.                                         *
!----------------------------------------------------------------------*
if (iOpt /= 0) then
  write(ErrMsg,*) 'Illegal option flag:',iOpt
  call SysAbendMsg('ffxRun',ErrMsg,' ')
end if
iRc = 0
!----------------------------------------------------------------------*
! Does the runfile exist? If not abend.                                *
!----------------------------------------------------------------------*
call f_inquire(RunName,ok)

! Do not return error for querying a runfile that does not exist,
! but rather return "not found", patch 6.7.263

!if (.not. ok) call SysAbendMsg('ffxRun','RunFile does not exist',' ')
if (.not. ok) then
  !write(6,*) ' Warning! In ffxRun: runfile does not exist!'
  iRc = rcNotFound
  nData = 0
  RecTyp = TypUnk
  return
end if
!----------------------------------------------------------------------*
! Open runfile.                                                        *
!----------------------------------------------------------------------*
call OpnRun(iRc,Lu,iOpt)
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
! Locate record.                                                       *
!----------------------------------------------------------------------*
item = -1
do i=1,nToc
  CmpLab1 = TocLab(i)
  CmpLab2 = Label
  call UpCase(CmpLab1)
  call UpCase(CmpLab2)
  if (CmpLab1 == CmpLab2) item = i
end do
if (item == -1) then
  iRc = rcNotFound
  nData = 0
  RecTyp = TypUnk
  call DaClos(Lu)
  return
end if
nData = TocLen(item)
RecTyp = TocTyp(item)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
call DaClos(Lu)

return

end subroutine ffxRun
