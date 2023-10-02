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

use casvb_global, only: release
use Constants, only: Zero
use Definitions, only: iwp

implicit none
#include "main_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: nvb_alloc
integer(kind=iwp), external :: mstackr_cvb, nvb_cvb

if (release(2)) call mfreer_cvb(lv(1))
release(2) = .true.
release(3) = .false.

! Note zeroing of ORBS and CVB:
lv(1) = mstackr_cvb(norb*norb)
work(lv(1):lv(1)+norb*norb-1) = Zero
! (MXNVB should be upper bound on KBASIS & KBASISCVB):
nvb_alloc = max(nvb_cvb(kbasiscvb),nvb_cvb(kbasis),mxnvb)
lv(2) = mstackr_cvb(nvb_alloc)
work(lv(2):lv(2)+nvb_alloc-1) = Zero

return

end subroutine chop2_cvb
