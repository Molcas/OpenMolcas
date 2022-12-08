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

use RunFile_data, only: Arr2RunHdr, icRd, IDRun, nHdrSz, NulPtr, RunHdr, RunName, VNRun
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: iRc, Lu
integer(kind=iwp), intent(in) :: iOpt
integer(kind=iwp) :: Arr(nHdrSz), iDisk
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
if (.not. ok) call SysAbendmsg('gxRdRun','RunFile does not exist',' ')
!----------------------------------------------------------------------*
! Open runfile and check that file is ok.                              *
!----------------------------------------------------------------------*
Lu = isFreeUnit(11)

RunHdr%ID = NulPtr
RunHdr%Ver = NulPtr
call DaName(Lu,RunName)
iDisk = 0
call iDaFile(Lu,icRd,Arr,nHdrSz,iDisk)
call Arr2RunHdr(Arr)
if (RunHdr%ID /= IDrun) then
  call DaClos(Lu)
  call SysFilemsg('gxWrRun','Wrong file type, not a RunFile',Lu,' ')
  call Abend()
end if
if (RunHdr%Ver /= VNrun) then
  call DaClos(Lu)
  call SysFilemsg('gxWrRun','Wrong version of RunFile',Lu,' ')
  call Abend()
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine OpnRun
