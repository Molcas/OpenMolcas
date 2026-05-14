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
! Copyright (C) 2024, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine wrToc(Lu)

use RunFile_data, only: icWr, lw, nToc, RunHdr, Toc
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: Lu
integer(kind=iwp) :: iDisk
integer(kind=iwp), allocatable :: Tmp(:)
character(len=lw), allocatable :: TmpLab(:)

! Prevent the compiler from creating temporary arrays
! (this means we have to create them ourselves, though)

call mma_allocate(Tmp,nToc,label='Tmp')
call mma_allocate(TmpLab,nToc,label='TmpLab')

TmpLab(:) = Toc(:)%Lab
iDisk = RunHdr%DaLab
call cDaFile(Lu,icWr,TmpLab,lw*nToc,iDisk)

Tmp(:) = Toc(:)%Ptr
iDisk = RunHdr%DaPtr
call iDaFile(Lu,icWr,Tmp,nToc,iDisk)

Tmp(:) = Toc(:)%Len
iDisk = RunHdr%DaLen
call iDaFile(Lu,icWr,Tmp,nToc,iDisk)

Tmp(:) = Toc(:)%MaxLen
iDisk = RunHdr%DaMaxLen
call iDaFile(Lu,icWr,Tmp,nToc,iDisk)

Tmp(:) = Toc(:)%Typ
iDisk = RunHdr%DaTyp
call iDaFile(Lu,icWr,Tmp,nToc,iDisk)

call mma_deallocate(Tmp)
call mma_deallocate(TmpLab)

end subroutine wrToc
