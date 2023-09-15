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

subroutine change6_cvb()

implicit real*8(a-h,o-z)
logical changed
! ... Change of dimensioning variables ...
logical, external :: chpcmp_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "rls_cvb.fh"
#include "davtune_cvb.fh"
#include "change6.fh"

changed = .false.
if (changed) call touch_cvb('CIFREE')

nprorb = norb*(norb-1)
if (strucopt) then
  nprvb = nvb
  npr = nprorb+nvb
else if (.not. ((imethod == 4) .or. (imethod == 6))) then
  nprvb = 0
  npr = nprorb
else
  npr = 0
end if
if (chpcmp_cvb(npr)) changed = .true.
if ((.not. ((imethod == 4) .or. (imethod == 6))) .and. (ifinish == 0)) then
  ! Standard 2nd-order procedure:
  icase = 1
else if ((imethod == 4) .and. (icrit == 1) .and. (ifinish == 0)) then
  ! Overlap-based Davidson
  icase = 2
else if ((imethod == 4) .and. (icrit == 2) .and. (ifinish == 0)) then
  ! Energy-based Davidson
  icase = 3
else if ((imethod == 6) .or. (ifinish == 1) .or. (ifinish == 2)) then
  ! No arrays needed
  icase = 4
else
  icase = 5
end if
if (chpcmp_cvb(icase)) changed = .true.
if (changed) call touch_cvb('MEM6')

return

end subroutine change6_cvb
