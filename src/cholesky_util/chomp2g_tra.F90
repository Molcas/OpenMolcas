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
! Copyright (C) 2004, Jonas Bostrom                                    *
!***********************************************************************

subroutine ChoMP2g_Tra(COrb1,COrb2,Diag,DoDiag,iMoType1,iMoType2)
!
! Jonas Bostrom, Dec. 2004.
!
! Purpose: transform Cholesky vectors to (pq) MO basis.

use Cholesky, only: nSym
use ChoMP2, only: nMoMo, nMoType
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: COrb1(*), COrb2(*)
real(kind=wp), intent(_OUT_) :: Diag(*)
logical(kind=iwp), intent(in) :: DoDiag
integer(kind=iwp) :: iMoType1, iMoType2
integer(kind=iwp) :: iSym, iVecType, kOffD, lw
real(kind=wp), allocatable :: TraMax(:)

! Allocate remaining memory.
! --------------------------

! Check what type of Cholesky vector to make (fro-occ, occ-occ.....)
iVecType = iMoType2+(iMoType1-1)*nMoType

call mma_maxDBLE(lw)
call mma_allocate(TraMax,lW,Label='TraMax')

kOffD = 1
do iSym=1,nSym

  ! Open files for MO vectors.
  ! --------------------------

  call ChoMP2_OpenF(1,1,iSym)

  ! Transform vectors.
  ! ------------------

  call ChoMP2g_Tra_1(COrb1,COrb2,Diag(kOffD),DoDiag,TraMax,lW,iSym,iMoType1,iMoType2)
  kOffD = kOffD+nMoMo(iSym,iVecType)

  ! Close files for MO vectors.
  ! ---------------------------

  call ChoMP2_OpenF(2,1,iSym)

end do

! Free memory.
! ------------

call mma_deallocate(TraMax)

end subroutine ChoMP2g_Tra
