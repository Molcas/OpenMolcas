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

subroutine change7_cvb()

use casvb_global, only: icase7, icrit, ifinish, imethod
use Definitions, only: iwp

implicit none
logical(kind=iwp) :: changed
logical(kind=iwp), external :: chpcmp_cvb ! ... Change of dimensioning variables ...

! Inconsequential work arrays
changed = .false.
if (((imethod /= 4) .and. (ifinish == 0)) .or. ((ifinish == 1) .or. (ifinish == 2))) then
  icase7 = 1
else if ((imethod == 4) .and. (icrit == 1) .and. (ifinish == 0)) then
  icase7 = 2
else if ((imethod == 4) .and. (icrit == 2) .and. (ifinish == 0)) then
  icase7 = 3
else
  icase7 = 4
end if
if (chpcmp_cvb(icase7)) changed = .true.
if (changed) call touch_cvb('MEM7')

return

end subroutine change7_cvb
