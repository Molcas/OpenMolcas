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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_Tra(COcc,CVir,Diag,DoDiag)
!
! Thomas Bondo Pedersen, Dec. 2004.
!
! Purpose: transform Cholesky vectors to (ai) MO basis.

use Cholesky, only: nSym
use ChoMP2, only: nT1am
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: COcc(*), CVir(*)
real(kind=wp), intent(inout) :: Diag(*)
logical(kind=iwp), intent(in) :: DoDiag
integer(kind=iwp) :: iSym, kOffD, lW
real(kind=wp), allocatable :: TraMax(:)

! Allocate remaining memory.
! --------------------------

call mma_maxDBLE(lw)
call mma_allocate(TraMax,lw,Label='TraMax')

kOffD = 1
do iSym=1,nSym

  ! Open files for MO vectors.
  ! --------------------------

  call ChoMP2_OpenF(1,1,iSym)

  ! Transform vectors.
  ! ------------------

  call ChoMP2_Tra_1(COcc,CVir,Diag(kOffD),DoDiag,TraMax,lW,iSym)
  if (DoDiag) kOffD = kOffD+nT1am(iSym)

  ! Close files for MO vectors.
  ! ---------------------------

  call ChoMP2_OpenF(2,1,iSym)

end do

! Free memory.
! ------------

call mma_deallocate(TraMax)

end subroutine ChoMP2_Tra
