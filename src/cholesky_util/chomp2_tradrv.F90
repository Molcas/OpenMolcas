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

subroutine ChoMP2_TraDrv(irc,CMO,Diag,DoDiag)
!
! Thomas Bondo Pedersen, Dec. 2004.
!
! Purpose: AO-to-MO (ai) transformation of Cholesky vectors
!          performed directly in reduced sets. This assumes
!          that the MP2 program has been appropriately initialized.

use ChoMP2, only: nAOVir, nT1AOT
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(in) :: CMO(*)
real(kind=wp), intent(inout) :: Diag(*)
logical(kind=iwp), intent(in) :: DoDiag
real(kind=wp), allocatable :: COcc(:), CVir(:)

irc = 0

! Reorder MO coefficients.
! ------------------------

call mma_allocate(COcc,nT1AOT(1),Label='COcc')
call mma_allocate(CVir,nAOVir(1),Label='CVir')
call ChoMP2_MOReOrd(CMO,COcc,CVir)

! Transform vectors.
! ------------------

call ChoMP2_Tra(COcc,CVir,Diag,DoDiag)

! Deallocate reordered MO coefficients.
! -------------------------------------

call mma_deallocate(CVir)
call mma_deallocate(COcc)

end subroutine ChoMP2_TraDrv
