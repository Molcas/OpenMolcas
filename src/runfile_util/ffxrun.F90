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

use RunFile_data, only: icRd, lw, nToc, rcNotFound, rcOK, RunHdr, RunName, Toc, TypUnk
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: iRc, nData, RecTyp
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(in) :: iOpt
integer(kind=iwp) :: i, iDisk, item, Lu
logical(kind=iwp) :: ok
character(len=lw) :: CmpLab1, CmpLab2
character(len=64) :: ErrMsg

!----------------------------------------------------------------------*
! Check that arguments are ok.                                         *
!----------------------------------------------------------------------*
if (iOpt /= 0) then
  write(ErrMsg,*) 'Illegal option flag:',iOpt
  call SysAbendMsg('ffxRun',ErrMsg,' ')
end if
iRc = rcOK
!----------------------------------------------------------------------*
! Does the runfile exist? If not abend.                                *
!----------------------------------------------------------------------*
call f_inquire(RunName,ok)

! Do not return error for querying a runfile that does not exist,
! but rather return "not found", patch 6.7.263

!if (.not. ok) call SysAbendMsg('ffxRun','RunFile does not exist',' ')
if (.not. ok) then
  !write(u6,*) ' Warning! In ffxRun: runfile does not exist!'
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
iDisk = RunHdr%DaLab
call cDaFile(Lu,icRd,Toc(:)%Lab,lw*nToc,iDisk)
iDisk = RunHdr%DaPtr
call iDaFile(Lu,icRd,Toc(:)%Ptr,nToc,iDisk)
iDisk = RunHdr%DaLen
call iDaFile(Lu,icRd,Toc(:)%Len,nToc,iDisk)
iDisk = RunHdr%DaMaxLen
call iDaFile(Lu,icRd,Toc(:)%MaxLen,nToc,iDisk)
iDisk = RunHdr%DaTyp
call iDaFile(Lu,icRd,Toc(:)%Typ,nToc,iDisk)
!----------------------------------------------------------------------*
! Locate record.                                                       *
!----------------------------------------------------------------------*
item = -1
do i=1,nToc
  CmpLab1 = Toc(i)%Lab
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
nData = Toc(item)%Len
RecTyp = Toc(item)%Typ
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
call DaClos(Lu)

return

end subroutine ffxRun
