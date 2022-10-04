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
! This routine prints the contents of the runfile.                     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************

subroutine DumpRun(iRc,iOpt)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iRc, iOpt
#include "runinfo.fh"
#include "runtypes.fh"
integer(kind=iwp) :: i, iDisk, Lu
character(len=64) :: ErrMsg

!----------------------------------------------------------------------*
! Check that arguments are ok.                                         *
!----------------------------------------------------------------------*
if (iOpt /= 0) then
  write(ErrMsg,*) 'Illegal option flag:',iOpt
  call SysAbendMsg('DumpRun',ErrMsg,' ')
end if
iRc = 0
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
! Print record information.                                            *
!----------------------------------------------------------------------*
write(u6,*)
write(u6,'(2a)') '------------------------------------------------------'
write(u6,'(a)') 'Contents in RunFile'
write(u6,'(2a)') '------------------------------------------------------'
write(u6,'(2a)') '  Slot        Label       Disk loc.   Field len.  Type'
write(u6,'(2a)') '  ----  ----------------  ----------  ----------  ----'
do i=1,nToc
  if (TocPtr(i) /= NulPtr) write(u6,'(i6,2x,a16,i12,2i12,i6)') i,TocLab(i),TocPtr(i),TocLen(i),TocMaxLen(i),TocTyp(i)
end do
write(u6,'(2a)') '------------------------------------------------------'
write(u6,*)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
call DaClos(Lu)

return

end subroutine DumpRun
