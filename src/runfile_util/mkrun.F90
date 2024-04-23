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

use RunFile_data, only: icWr, IDRun, lw, nHdrSz, nToc, NulPtr, RunHdr, RunHdr2Arr, RunName, Toc, TypUnk, VNRun
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp
!use Definitions, only: u6
use Para_Info, only: nProcs

implicit none
integer(kind=iwp), intent(out) :: iRc
integer(kind=iwp), intent(in) :: iOpt
integer(kind=iwp) :: Hdr(nHdrSz), iAllow, iDisk, Lu
logical(kind=iwp) :: ok
character(len=64) :: ErrMsg
integer(kind=iwp), allocatable :: Tmp(:)
character(len=lw), allocatable :: TmpLab(:)
integer(kind=iwp), external :: isFreeUnit

!----------------------------------------------------------------------*
! Check that arguments are ok.                                         *
!----------------------------------------------------------------------*
iAllow = not(ibset(0,0))
if (iand(iOpt,iAllow) /= 0) then
  write(ErrMsg,*) 'Illegal option flag:',iOpt
  call SysAbendMsg('MkRun',ErrMsg,' ')
end if
iRc = 0
!----------------------------------------------------------------------*
! Optionally do not create.                                            *
!----------------------------------------------------------------------*
if (btest(iOpt,0)) then
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

RunHdr%ID = IDrun
RunHdr%Ver = VNrun
RunHdr%Next = 0
RunHdr%Items = 0
RunHdr%nProcs=nProcs
call DaName(Lu,RunName)
iDisk = 0
call RunHdr2Arr(Hdr)
call iDaFile(Lu,icWr,Hdr,nHdrSz,iDisk)
RunHdr%Next = iDisk
iDisk = 0
call RunHdr2Arr(Hdr)
call iDaFile(Lu,icWr,Hdr,nHdrSz,iDisk)

iDisk = RunHdr%Next
call mma_allocate(Tmp,nToc,label='Tmp')
call mma_allocate(TmpLab,nToc,label='TmpLab')

TmpLab(:) = 'Empty'
RunHdr%DaLab = iDisk
call cDaFile(Lu,icWr,TmpLab,lw*nToc,iDisk)
Toc(:)%Lab = TmpLab(:)

Tmp(:) = NulPtr
RunHdr%DaPtr = iDisk
call iDaFile(Lu,icWr,Tmp,nToc,iDisk)
Toc(:)%Ptr = Tmp(:)

Tmp(:) = 0
RunHdr%DaLen = iDisk
call iDaFile(Lu,icWr,Tmp,nToc,iDisk)
Toc(:)%Len = Tmp(:)

RunHdr%DaMaxLen = iDisk
call iDaFile(Lu,icWr,Tmp,nToc,iDisk)
Toc(:)%MaxLen = Tmp(:)

Tmp(:) = TypUnk
RunHdr%DaTyp = iDisk
call iDaFile(Lu,icWr,Tmp,nToc,iDisk)
Toc(:)%Typ = Tmp(:)

call mma_deallocate(Tmp)
call mma_deallocate(TmpLab)
RunHdr%Next = iDisk
iDisk = 0
call RunHdr2Arr(Hdr)
call iDaFile(Lu,icWr,Hdr,nHdrSz,iDisk)

call DaClos(Lu)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine MkRun
