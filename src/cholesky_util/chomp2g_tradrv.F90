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
! Copyright (C) 2010, Jonas Bostrom                                    *
!***********************************************************************

subroutine ChoMP2g_TraDrv(irc,CMO,Diag,DoDiag)
!
! Jonas Bostrom, Jan. 2010. (modified from ChoMP2_TraDrv)
!
! Purpose: AO-to-MO (pq) transformation of Cholesky vectors
!          performed directly in reduced sets. This assumes
!          that the MP2 program has been appropriately initialized.

use Cholesky, only: nSym
use ChoMP2, only: nAdrOff, nMoAo, nMoType
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(in) :: CMO(*)
real(kind=wp), intent(_OUT_) :: Diag(*)
logical(kind=iwp), intent(in) :: DoDiag
integer(kind=iwp) :: i, iMoType, iSym, jMOType, l_COrb, nProdType
logical(kind=iwp) :: DoDiag_
real(kind=wp), allocatable :: COrb1(:), COrb2(:)

#include "macros.fh"
unused_var(DoDiag)

irc = 0

! Reorder MO coefficients.
! ------------------------

DoDiag_ = .false.
nProdType = nMOType**2
l_COrb = 0
nAdrOff(1:nSym) = 0

do iSym=1,nSym
  do i=1,nProdType
    l_COrb = max(l_COrb,nMoAo(iSym,i))
  end do
end do
call mma_allocate(COrb1,l_COrb,Label='COrb1')
call mma_allocate(COrb2,l_COrb,Label='COrb2')

DoDiag_ = .true.
call ChoMP2g_MOReOrd(CMO,COrb1,COrb2,2,3)
call ChoMP2g_Tra(COrb1,COrb2,Diag,DoDiag_,2,3)
DoDiag_ = .false.
do iMoType=1,3
  do jMOType=1,3
    if ((iMoType == 2) .and. (jMoType == 3)) cycle

    call ChoMP2g_MOReOrd(CMO,COrb1,COrb2,iMOType,jMOType)
    ! Transform vectors.
    ! ------------------
    call ChoMP2g_Tra(COrb1,COrb2,Diag,DoDiag_,iMoType,jMoType)

  end do
end do

! Deallocate reordered MO coefficients.
! -------------------------------------

call mma_deallocate(COrb2)
call mma_deallocate(COrb1)

end subroutine ChoMP2g_TraDrv
