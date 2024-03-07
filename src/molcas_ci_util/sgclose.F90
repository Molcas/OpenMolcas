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
Subroutine SGClose(SGS)
use stdalloc, only: mma_deallocate
use Struct, only: SGStruct
Implicit None
Type (SGStruct) SGS

If (Allocated(SGS%ISM)) Call mma_deallocate(SGS%ISM)
If (Allocated(SGS%DRT0)) Call mma_deallocate(SGS%DRT0)
If (Allocated(SGS%DOWN0)) Call mma_deallocate(SGS%DOWN0)
If (Allocated(SGS%DRT)) Call mma_deallocate(SGS%DRT)
If (Allocated(SGS%DOWN)) Call mma_deallocate(SGS%DOWN)
If (Allocated(SGS%UP)) Call mma_deallocate(SGS%UP)
If (Allocated(SGS%MAW)) Call mma_deallocate(SGS%MAW)
If (Allocated(SGS%LTV)) Call mma_deallocate(SGS%LTV)
If (Allocated(SGS%DAW)) Call mma_deallocate(SGS%DAW)
If (Allocated(SGS%RAW)) Call mma_deallocate(SGS%RAW)
If (Allocated(SGS%SCR)) Call mma_deallocate(SGS%SCR)
If (Allocated(SGS%Ver)) Call mma_deallocate(SGS%Ver)
SGS%DRTP => Null()
SGS%DOWNP => Null()

end Subroutine SGClose
