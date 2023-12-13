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

use RunFile_data, only: icRd, lw, nToc, NulPtr, RunHdr, Toc
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iRc
integer(kind=iwp), intent(in) :: iOpt
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
! Print record information.                                            *
!----------------------------------------------------------------------*
write(u6,*)
write(u6,'(2a)') '------------------------------------------------------'
write(u6,'(a)') 'Contents in RunFile'
write(u6,'(2a)') '------------------------------------------------------'
write(u6,'(2a)') '  Slot        Label       Disk loc.   Field len.  Type'
write(u6,'(2a)') '  ----  ----------------  ----------  ----------  ----'
do i=1,nToc
  if (Toc(i)%Ptr /= NulPtr) write(u6,'(i6,2x,a,i12,2i12,i6)') i,Toc(i)%Lab,Toc(i)%Ptr,Toc(i)%Len,Toc(i)%MaxLen,Toc(i)%Typ
end do
write(u6,'(2a)') '------------------------------------------------------'
write(u6,*)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
call DaClos(Lu)

return

end subroutine DumpRun
