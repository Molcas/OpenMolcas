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
! This routine creates a runfile.                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! iOpt: Bitswitch                                                      *
!    1 -- Do not create if exist.                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************

subroutine MkRun(iRc,iOpt)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iRc, iOpt
#include "runinfo.fh"
#include "runtypes.fh"
integer(kind=iwp) :: i, iAllow, iDisk, Lu
logical(kind=iwp) :: ok
character(len=64) :: ErrMsg
integer(kind=iwp), external :: isFreeUnit

!----------------------------------------------------------------------*
! Check that arguments are ok.                                         *
!----------------------------------------------------------------------*
iAllow = -1
iAllow = ieor(iAllow,1)
if (iand(iOpt,iAllow) /= 0) then
  write(ErrMsg,*) 'Illegal option flag:',iOpt
  call SysAbendMsg('MkRun',ErrMsg,' ')
end if
iRc = 0
!----------------------------------------------------------------------*
! Optionally do not create.                                            *
!----------------------------------------------------------------------*
if (iand(iOpt,1) /= 0) then
  call f_inquire(RunName,ok)
  if (ok) then
    !write(u6,*) '*** NOT creating runfile ',RunName
    return
  end if
end if
!write(u6,*) '*** Creating runfile ',RunName
!----------------------------------------------------------------------*
! Create it.                                                           *
!----------------------------------------------------------------------*
Lu = 11
Lu = isFreeUnit(Lu)

RunHdr(ipID) = IDrun
RunHdr(ipVer) = VNrun
RunHdr(ipNext) = 0
RunHdr(ipItems) = 0
call DaName(Lu,RunName)
iDisk = 0
call iDaFile(Lu,icWr,RunHdr,nHdrSz,iDisk)
RunHdr(ipNext) = iDisk
iDisk = 0
call iDaFile(Lu,icWr,RunHdr,nHdrSz,iDisk)

iDisk = RunHdr(ipNext)
do i=1,nToc
  TocLab(i) = 'Empty   '
  TocPtr(i) = NulPtr
  TocLen(i) = 0
  TocMaxLen(i) = 0
  TocTyp(i) = TypUnk
end do
RunHdr(ipDaLab) = iDisk
call cDaFile(Lu,icWr,TocLab,16*nToc,iDisk)
RunHdr(ipDaPtr) = iDisk
call iDaFile(Lu,icWr,TocPtr,nToc,iDisk)
RunHdr(ipDaLen) = iDisk
call iDaFile(Lu,icWr,TocLen,nToc,iDisk)
RunHdr(ipDaMaxLen) = iDisk
call iDaFile(Lu,icWr,TocMaxLen,nToc,iDisk)
RunHdr(ipDaTyp) = iDisk
call iDaFile(Lu,icWr,TocTyp,nToc,iDisk)
RunHdr(ipNext) = iDisk
iDisk = 0
call iDaFile(Lu,icWr,RunHdr,nHdrSz,iDisk)

call DaClos(Lu)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine MkRun
