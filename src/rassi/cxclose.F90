!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
Subroutine CXClose(CIS,EXS)
use stdalloc, only: mma_deallocate
use Struct, only: CIStruct, EXStruct
Type (CIStruct) CIS
Type (ExStruct) ExS
IF (Allocated(EXS%MVL)) Call mma_deallocate(EXS%MVL)
IF (Allocated(EXS%MVR)) Call mma_deallocate(EXS%MVR)
IF (Allocated(CIS%NOCSF)) Call mma_deallocate(CIS%NOCSF)
IF (Allocated(CIS%IOCSF)) Call mma_deallocate(CIS%IOCSF)
IF (Allocated(CIS%NOW)) Call mma_deallocate(CIS%NOW)
IF (Allocated(CIS%IOW)) Call mma_deallocate(CIS%IOW)
IF (Allocated(EXS%NOCP)) Call mma_deallocate(EXS%NOCP)
IF (Allocated(EXS%IOCP)) Call mma_deallocate(EXS%IOCP)
IF (Allocated(CIS%NCSF)) Call mma_deallocate(CIS%NCSF)
IF (Allocated(CIS%ICase)) Call mma_deallocate(CIS%ICase)
IF (Allocated(EXS%ICoup)) Call mma_deallocate(EXS%ICoup)
IF (Allocated(EXS%VTab)) Call mma_deallocate(EXS%VTab)
IF (Allocated(EXS%SGTMP)) Call mma_deallocate(EXS%SGTMP)

end subroutine CXClose
