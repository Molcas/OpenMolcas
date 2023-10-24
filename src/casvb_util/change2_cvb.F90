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

subroutine change2_cvb()

use casvb_global, only: kbasis, kbasiscvb, mxnvb, norb
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nvb_alloc
logical(kind=iwp) :: changed
integer(kind=iwp), external :: nvb_cvb
logical(kind=iwp), external :: chpcmp_cvb ! ... Change of dimensioning variables ...

changed = .false.
! VB wavefunction
! (MXNVB should be upper bound on KBASIS & KBASISCVB):
nvb_alloc = max(nvb_cvb(kbasiscvb),nvb_cvb(kbasis),mxnvb)
if (chpcmp_cvb(norb)) changed = .true.
if (chpcmp_cvb(nvb_alloc)) changed = .true.
if (changed) call touch_cvb('MEM2')

return

end subroutine change2_cvb
