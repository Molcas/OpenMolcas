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
Call mma_deallocate(EXS%MVL)
Call mma_deallocate(EXS%MVR)
Call mma_deallocate(CIS%NOCSF)
Call mma_deallocate(CIS%IOCSF)
Call mma_deallocate(CIS%NOW)
Call mma_deallocate(CIS%IOW)
Call mma_deallocate(EXS%NOCP)
Call mma_deallocate(EXS%IOCP)
Call mma_deallocate(CIS%NCSF)
Call mma_deallocate(CIS%ICase)
Call mma_deallocate(EXS%ICoup)
Call mma_deallocate(EXS%VTab)

end subroutine CXClose
