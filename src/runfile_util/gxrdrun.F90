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
! This routine reads a record from the runfile.                        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************

subroutine gxRdRun(iRc,Label,cData,nData,iOpt,RecTyp)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iRc, nData, iOpt, RecTyp
character(len=*) :: Label
character :: cData(*)
#include "runinfo.fh"
#include "runtypes.fh"
integer(kind=iwp) :: DataAdr, i, iDisk, item, Lu
logical(kind=iwp) :: ok
character(len=64) :: ErrMsg
character(len=16) :: CmpLab1, CmpLab2

!----------------------------------------------------------------------*
! Check that arguments are ok.                                         *
!----------------------------------------------------------------------*
ok = .false.
if (RecTyp == TypInt) ok = .true.
if (RecTyp == TypDbl) ok = .true.
if (RecTyp == TypStr) ok = .true.
if (RecTyp == TypLgl) ok = .true.
if (.not. ok) call SysAbendMsg('gxRdRun','Argument RecTyp is of wrong type','Aborting')
if (nData < 0) call SysAbendMsg('gxRdRun','Number of data items less than zero','Aborting')
if (iOpt /= 0) then
  write(ErrMsg,*) 'Illegal option flag:',iOpt
  call SysAbendMsg('gxRdRun',ErrMsg,' ')
end if
iRc = 0
!----------------------------------------------------------------------*
! Does the runfile exist? If not abort.                                *
!----------------------------------------------------------------------*
call f_inquire(RunName,ok)
if (.not. ok) call SysFilemsg('gxRdRun','RunFile does not exist',Lu,' ')
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
! Find field.                                                          *
!----------------------------------------------------------------------*
item = -1
do i=1,nToc
  CmpLab1 = TocLab(i)
  CmpLab2 = Label
  !call Upcase(CmpLab1)
  !call Upcase(CmpLab2)
  if (CmpLab1 == CmpLab2) item = i
end do
if (item == -1) then
  call DaClos(Lu)
  write(ErrMsg,'(a,a)') 'Record not found in runfile: ',Label
  call SysFilemsg('gxRdRun',ErrMsg,Lu,' ')
end if
DataAdr = TocPtr(item)
!----------------------------------------------------------------------*
! Read data from runfile.                                              *
!----------------------------------------------------------------------*
iDisk = DataAdr
call gzRWRun(Lu,icRd,cData,nData,iDisk,RecTyp)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
call DaClos(Lu)

return

end subroutine gxRdRun
