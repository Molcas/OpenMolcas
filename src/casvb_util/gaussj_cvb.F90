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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

!***********************************************************************
!*                                                                     *
!*  GAUSSJ    := Define sequence of simple updates from orb transf.    *
!*                                                                     *
!***********************************************************************
subroutine gaussj_cvb(orbs,gjorb)

use casvb_global, only: gjorb_type, norb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: orbs(norb,norb)
type(gjorb_type), intent(_OUT_) :: gjorb
integer(kind=iwp) :: i
integer(kind=iwp), allocatable :: lrow(:)
real(kind=wp), allocatable :: a(:,:)

call mma_allocate(a,norb,norb,label='a')
call mma_allocate(lrow,norb,label='lrow')
a(:,:) = orbs(:,:)
call gaussj2_cvb(a,lrow,gjorb%i1,gjorb%i2,gjorb%r,norb)
lrow(:) = gjorb%i1(:)
do i=1,norb
  gjorb%i1(lrow(i)) = i
end do
call mma_deallocate(a)
call mma_deallocate(lrow)

return

end subroutine gaussj_cvb
