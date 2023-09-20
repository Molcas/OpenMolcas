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

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

if (release(2)) call mfreer_cvb(lv(1))
release(2) = .true.
release(3) = .false.

! Note zeroing of ORBS and CVB:
lv(1) = mstackrz_cvb(norb*norb)
! (MXNVB should be upper bound on KBASIS & KBASISCVB):
nvb_alloc = max(nvb_cvb(kbasiscvb),nvb_cvb(kbasis),mxnvb)
lv(2) = mstackrz_cvb(nvb_alloc)

return

end subroutine chop2_cvb
