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

use RunFile_data, only: icRd, lw, nToc, RunHdr, RunName, Toc, TypDbl, TypInt, TypLgl, TypStr
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: iRc
character(len=*), intent(in) :: Label
character, intent(_OUT_) :: cData(*)
integer(kind=iwp), intent(in) :: nData, iOpt, RecTyp
integer(kind=iwp) :: DataAdr, i, iDisk, item, Lu
logical(kind=iwp) :: ok
character(len=lw) :: CmpLab1, CmpLab2
character(len=64) :: ErrMsg

!----------------------------------------------------------------------*
! Check that arguments are ok.                                         *
!----------------------------------------------------------------------*
select case (RecTyp)
  case (TypInt,TypDbl,TypStr,TypLgl)
    !continue ! ok
  case default
    call SysAbendMsg('gxRdRun','Argument RecTyp is of wrong type','Aborting')
end select
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
if (.not. ok) call SysAbendmsg('gxRdRun','RunFile does not exist',' ')
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
! Find field.                                                          *
!----------------------------------------------------------------------*
item = -1
do i=1,nToc
  CmpLab1 = Toc(i)%Lab
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
DataAdr = Toc(item)%Ptr
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
