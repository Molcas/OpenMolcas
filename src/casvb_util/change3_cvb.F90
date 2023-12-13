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

subroutine change3_cvb()

use casvb_global, only: kbasis, kbasiscvb
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: kmost
logical(kind=iwp) :: changed
logical(kind=iwp), external :: chpcmp_cvb ! ... Change of dimensioning variables ...

changed = .false.
! Spin functions coefficients (BIKCOF) + inverse (AIKCOF)
!  Get KBASISCVB if we don't know it already (eqv. to GETGUESS later):
!  Need to reserve enough space for both KBASIS and KBASISCVB
!  --> figure out which one needs most:
if (((kbasis > 2) .and. (kbasis /= 6)) .or. ((kbasiscvb > 2) .and. (kbasiscvb /= 6))) then
  kmost = 3
else if ((kbasis <= 2) .or. (kbasiscvb <= 2)) then
  kmost = 1
else
  kmost = 6
end if
if (chpcmp_cvb(kmost)) changed = .true.
if (changed) call touch_cvb('MEM3')

return

end subroutine change3_cvb
