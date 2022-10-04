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
! This routine opens the runfile.                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************

subroutine OpnRun(iRc,Lu,iOpt)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iRc, Lu, iOpt
#include "runinfo.fh"
integer(kind=iwp) :: iDisk
logical(kind=iwp) :: ok
character(len=64) :: ErrMsg
integer(kind=iwp), external :: isFreeUnit

!----------------------------------------------------------------------*
! Check that arguments are ok.                                         *
!----------------------------------------------------------------------*
if (iOpt /= 0) then
  write(ErrMsg,*) 'Illegal option flag:',iOpt
  call SysAbendMsg('OpnRun',ErrMsg,' ')
end if
iRc = 0
!----------------------------------------------------------------------*
! Does the runfile exist? If not abort.                                *
!----------------------------------------------------------------------*
call f_inquire(RunName,ok)
if (.not. ok) call SysFilemsg('gxRdRun','RunFile does not exist',Lu,' ')
!----------------------------------------------------------------------*
! Open runfile and check that file is ok.                              *
!----------------------------------------------------------------------*
Lu = 11
Lu = isFreeUnit(Lu)

RunHdr(ipID) = -77
RunHdr(ipVer) = -77
call DaName(Lu,RunName)
iDisk = 0
call iDaFile(Lu,icRd,RunHdr,nHdrSz,iDisk)
if (RunHdr(ipID) /= IDrun) then
  call DaClos(Lu)
  call SysFilemsg('gxWrRun','Wrong file type, not a RunFile',Lu,' ')
  call Abend()
end if
if (RunHdr(ipVer) /= VNrun) then
  call DaClos(Lu)
  call SysFilemsg('gxWrRun','Wrong version of RunFile',Lu,' ')
  call Abend()
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine OpnRun
