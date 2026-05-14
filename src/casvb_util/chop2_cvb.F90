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

subroutine chop2_cvb()

use casvb_global, only: cvb, kbasis, kbasiscvb, mxnvb, norb, orbs, release
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nvb_alloc
integer(kind=iwp), external :: nvb_cvb

if (release(2)) then
  call mma_deallocate(orbs)
  call mma_deallocate(cvb)
end if
release(2) = .true.
release(3) = .false.

! Note zeroing of ORBS and CVB:
call mma_allocate(orbs,norb,norb,label='orbs')
orbs(:,:) = Zero
! (MXNVB should be upper bound on KBASIS & KBASISCVB):
nvb_alloc = max(nvb_cvb(kbasiscvb),nvb_cvb(kbasis),mxnvb)
call mma_allocate(cvb,nvb_alloc,label='cvb')
cvb(:) = Zero

return

end subroutine chop2_cvb
