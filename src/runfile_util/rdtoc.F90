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

subroutine rdToc(Lu)

use RunFile_data, only: icRd, lw, nToc, RunHdr, Toc
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

iDisk = RunHdr%DaLab
call cDaFile(Lu,icRd,TmpLab,lw*nToc,iDisk)
Toc(:)%Lab = TmpLab(:)

iDisk = RunHdr%DaPtr
call iDaFile(Lu,icRd,Tmp,nToc,iDisk)
Toc(:)%Ptr = Tmp(:)

iDisk = RunHdr%DaLen
call iDaFile(Lu,icRd,Tmp,nToc,iDisk)
Toc(:)%Len = Tmp(:)

iDisk = RunHdr%DaMaxLen
call iDaFile(Lu,icRd,Tmp,nToc,iDisk)
Toc(:)%MaxLen = Tmp(:)

iDisk = RunHdr%DaTyp
call iDaFile(Lu,icRd,Tmp,nToc,iDisk)
Toc(:)%Typ = Tmp(:)

call mma_deallocate(Tmp)
call mma_deallocate(TmpLab)

end subroutine rdToc
