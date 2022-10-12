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

subroutine gxWrRun(iRc,Label,cData,nData,iOpt,RecTyp)

use RunFile_data, only: icRd, icWr, lw, nHdrSz, nToc, NulPtr, RunHdr, RunHdr2Arr, RunName, Toc, TypDbl, TypInt, TypLgl, TypStr, &
                        TypUnk
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: iRc
character(len=*), intent(in) :: Label
character, intent(in) :: cData(*)
integer(kind=iwp), intent(in) :: nData, iOpt, RecTyp
integer(kind=iwp) :: Hdr(nHdrSz), DataAdr, i, iDisk, item, Lu, NewLen
logical(kind=iwp) :: ok, remove
character(len=64) :: ErrMsg

!----------------------------------------------------------------------*
! Check that arguments are ok.                                         *
!----------------------------------------------------------------------*
DataAdr = -99999
select case (RecTyp)
  case (TypInt,TypDbl,TypStr,TypLgl)
    !continue ! ok
  case default
    call SysAbendMsg('gxWrRun','Argument RecTyp is of wrong type','Aborting')
end select
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
if (RunHdr%Items >= nToc) then
  call DaClos(Lu)
  call SysFilemsg('gxWrRun','Ran out of ToC record in RunFile',Lu,' ')
  call Abend()
end if
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
! Reuse old field?                                                     *
!----------------------------------------------------------------------*
item = -1
do i=1,nToc
  if (Toc(i)%Lab == Label) item = i
end do
NewLen = 0
if (item /= -1) then
  remove = .false.
  if (Toc(item)%Typ /= RecTyp) remove = .true.
  if (Toc(item)%MaxLen < nData) remove = .true.
  if (remove) then
    !write (u6,*) '*******************************************'
    !write (u6,'(a,a,a,i10)') 'Label=',Label,' expands in RUNFILE with size=',nData
    !write (u6,*) '*******************************************'
    !call Abend()
    Toc(item)%Lab = 'Empty   '
    Toc(item)%Ptr = NulPtr
    Toc(item)%Len = 0
    Toc(item)%Typ = TypUnk
    item = -1
  else
    DataAdr = Toc(item)%Ptr
    NewLen = Toc(item)%Len
  end if
  RunHdr%Items = RunHdr%Items-1
end if
!----------------------------------------------------------------------*
! Use new field?                                                       *
!----------------------------------------------------------------------*
if (item == -1) then
  do i=nToc,1,-1
    if (Toc(i)%Ptr == NulPtr) item = i
  end do
  if (item == -1) then
    call DaClos(Lu)
    call SysFilemsg('gxWrRun','Internal inconsistency handling RunFile',Lu,' ')
    call Abend()
  end if
  DataAdr = RunHdr%Next
end if
!----------------------------------------------------------------------*
! Write data to runfile and update header.                             *
!----------------------------------------------------------------------*
RunHdr%Items = RunHdr%Items+1
Toc(item)%Lab = Label
Toc(item)%Ptr = DataAdr
Toc(item)%Len = nData
Toc(item)%MaxLen = max(NewLen,nData)
Toc(item)%Typ = RecTyp

iDisk = DataAdr
call gzRWRun(Lu,icWr,cData,nData,iDisk,RecTyp)

if (iDisk > RunHdr%Next) RunHdr%Next = iDisk
iDisk = 0
call RunHdr2Arr(Hdr)
call iDaFile(Lu,icWr,Hdr,nHdrSz,iDisk)
iDisk = RunHdr%DaLab
call cDaFile(Lu,icWr,Toc(:)%Lab,lw*nToc,iDisk)
iDisk = RunHdr%DaPtr
call iDaFile(Lu,icWr,Toc(:)%Ptr,nToc,iDisk)
iDisk = RunHdr%DaLen
call iDaFile(Lu,icWr,Toc(:)%Len,nToc,iDisk)
iDisk = RunHdr%DaMaxLen
call iDaFile(Lu,icWr,Toc(:)%MaxLen,nToc,iDisk)
iDisk = RunHdr%DaTyp
call iDaFile(Lu,icWr,Toc(:)%Typ,nToc,iDisk)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
call DaClos(Lu)

return

end subroutine gxWrRun
