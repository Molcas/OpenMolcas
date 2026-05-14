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

subroutine makegjorbs_cvb(orbs)
! Construct Gauss-Jordan factorizations of ORBS, ORBS transpose,
! and overlap matrix corresonding to ORBS:

use casvb_global, only: gjorb, gjorb2, gjorb3, norb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp

implicit none
real(kind=wp), intent(in) :: orbs(norb,norb)
real(kind=wp), allocatable :: owrk(:,:)

call mma_allocate(owrk,norb,norb,label='owrk')

call gaussj_cvb(orbs,gjorb)

call trnsps(norb,norb,orbs,owrk)
call gaussj_cvb(owrk,gjorb2)

call mxattb_cvb(orbs,orbs,norb,norb,norb,owrk)
call gaussj_cvb(owrk,gjorb3)

call mma_deallocate(owrk)

return

end subroutine makegjorbs_cvb
