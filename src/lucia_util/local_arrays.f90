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
Module Local_Arrays
use stdalloc, only: mma_allocate, mma_deallocate

Private

! CLBT : Length of each Batch ( in blocks)
! CLEBT : Length of each Batch ( in elements)
! CI1BT : Length of each block
! CIBT  : Info on each block
! CBLTP : BLock type for each symmetry

Integer, Allocatable, Public:: CLBT(:), CLEBT(:), CI1BT(:), CIBT(:)
Integer, Allocatable, Public, Target:: CBLTP(:)

Public :: Allocate_Local_Arrays, Deallocate_Local_Arrays
Contains

Subroutine Allocate_Local_Arrays(MXNTTS,NSMST)
Integer :: MXNTTS,NSMST
If (Allocated(CLBT)) Call Abend()
CALL mma_allocate(CLBT ,MXNTTS,Label='CLBT',safe='*')
CALL mma_allocate(CLEBT,MXNTTS,Label='CLEBT',safe='*')
CALL mma_allocate(CI1BT,MXNTTS,Label='CI1BT',safe='*')
CALL mma_allocate(CIBT ,8*MXNTTS,Label='CIBT',safe='*')
CALL mma_allocate(CBLTP,NSMST,Label='CBLTP',safe='*')
End Subroutine Allocate_Local_Arrays

Subroutine Deallocate_Local_Arrays()
CALL mma_deallocate(CLBT,safe='*')
CALL mma_deallocate(CLEBT,safe='*')
CALL mma_deallocate(CI1BT,safe='*')
CALL mma_deallocate(CIBT,safe='*')
CALL mma_deallocate(CBLTP,safe='*')
End Subroutine Deallocate_Local_Arrays
End Module Local_Arrays
