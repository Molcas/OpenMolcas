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
! This program converts a runfile into an ascii file.                  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
! Written: July 2003                                                   *
!          Lund University, Sweden                                     *
!                                                                      *
!***********************************************************************

program RF2Asc

use RunFile_data, only: lw, nToc, NulPtr, RunHdr, Toc, TypDbl, TypInt, TypLgl, TypStr
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, iDisk, iOpt, iRc, Lu, RUNASCII
logical(kind=iwp) :: Found
integer(kind=iwp), allocatable :: iBuf(:)
real(kind=wp), allocatable :: dBuf(:)
character, allocatable :: cBuf(:)
integer(kind=iwp), external :: isFreeUnit

call IniMem()
call Init_LinAlg()
call PrgmInit('RF2Asc')
!----------------------------------------------------------------------*
! Open runfile and check that file is ok.                              *
!----------------------------------------------------------------------*
call f_inquire('RUNFILE',Found)
if (.not. Found) then
  call WarningMessage(2,'RUNFILE not found')
  stop
end if
call NameRun('RUNFILE')
iOpt = 0
iRc = 0
call OpnRun(iRc,Lu,iOpt)
!----------------------------------------------------------------------*
! Read the ToC                                                         *
!----------------------------------------------------------------------*
iDisk = RunHdr%DaLab
call cDaFile(Lu,2,Toc(:)%Lab,lw*nToc,iDisk)
iDisk = RunHdr%DaPtr
call iDaFile(Lu,2,Toc(:)%Ptr,nToc,iDisk)
iDisk = RunHdr%DaLen
call iDaFile(Lu,2,Toc(:)%Len,nToc,iDisk)
iDisk = RunHdr%DaTyp
call iDaFile(Lu,2,Toc(:)%Typ,nToc,iDisk)
!----------------------------------------------------------------------*
! Open output file.                                                    *
!----------------------------------------------------------------------*
RUNASCII = isFreeUnit(15)
call MOLCAS_OPEN(RUNASCII,'RUNASCII')
!----------------------------------------------------------------------*
! Print header information.                                            *
!----------------------------------------------------------------------*
write(RUNASCII,'(a)') '# Header information'
write(RUNASCII,'(a,i15)') '# ID                      ',RunHdr%ID
write(RUNASCII,'(a,i15)') '# Version                 ',RunHdr%Ver
write(RUNASCII,'(a,i15)') '# Next free address       ',RunHdr%Next
write(RUNASCII,'(a,i15)') '# Number of items         ',RunHdr%Items
write(RUNASCII,'(a,i15)') '# Address to ToC labels   ',RunHdr%DaLab
write(RUNASCII,'(a,i15)') '# Address to ToC pointers ',RunHdr%DaPtr
write(RUNASCII,'(a,i15)') '# Address to ToC lengths  ',RunHdr%DaLen
write(RUNASCII,'(a,i15)') '# Address to ToC types    ',RunHdr%DaTyp
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
do i=1,nToc
  if (Toc(i)%Ptr == NulPtr) cycle
  write(RUNASCII,'(3a)') '<',Toc(i)%Lab,'>'
  write(RUNASCII,*) Toc(i)%Len,Toc(i)%Typ
  select case (Toc(i)%Typ)
    case (TypDbl)
      iDisk = Toc(i)%Ptr
      call mma_allocate(dBuf,Toc(i)%Len,Label='dBuf')
      call dDaFile(Lu,2,dBuf,Toc(i)%Len,iDisk)
      write(RUNASCII,'(4d26.18)') dBuf(1:Toc(i)%Len)
      call mma_deallocate(dBuf)
    case (TypInt)
      iDisk = Toc(i)%Ptr
      call mma_allocate(iBuf,Toc(i)%Len,Label='iBuf')
      call iDaFile(Lu,2,iBuf,Toc(i)%Len,iDisk)
      write(RUNASCII,*) iBuf(1:Toc(i)%Len)
      call mma_deallocate(iBuf)
    case (TypStr)
      iDisk = Toc(i)%Ptr
      call mma_allocate(cBuf,Toc(i)%Len,Label='cBuf')
      call cDaFile(Lu,2,cBuf,Toc(i)%Len,iDisk)
      write(RUNASCII,'(64a1)') cBuf(1:Toc(i)%Len)
      call mma_deallocate(cBuf)
    case (TypLgl)
      write(u6,*) 'Cannot handle type logical'
      exit
  end select
end do
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
close(RUNASCII)

call GetMem('Finish','LIST','REAL',I,0)
call GetMem('Finish','TERM','REAL',I,0)

end program RF2Asc
