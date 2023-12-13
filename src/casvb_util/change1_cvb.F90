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

subroutine change1_cvb()

use casvb_global, only: kbasis, mnion, mxion, nalf, nbet, nconf, ndetvb, nel, norb, nvb
use Definitions, only: iwp

implicit none
logical(kind=iwp) :: changed
integer(kind=iwp), external :: nvb_cvb
logical(kind=iwp), external :: chpcmp_cvb ! ... Change of dimensioning variables ...

! Arrays for determinant handling (E_ij) and definition of VB wavefunction:
changed = .false.
if (chpcmp_cvb(norb)) changed = .true.
if (chpcmp_cvb(nalf)) changed = .true.
if (chpcmp_cvb(nbet)) changed = .true.
if (chpcmp_cvb(nel)) changed = .true.

if (changed) call touch_cvb('CASPRINT')

if (chpcmp_cvb(nconf)) changed = .true.

if (.not. changed) call cnfchk_cvb()
nvb = nvb_cvb(kbasis)

if (chpcmp_cvb(ndetvb)) changed = .true.
if (chpcmp_cvb(mxion)) changed = .true.
if (chpcmp_cvb(mnion)) changed = .true.
if (changed) call touch_cvb('MEM1')

return

end subroutine change1_cvb
