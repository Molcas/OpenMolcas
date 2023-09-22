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

subroutine stat1_cvb()

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
#include "print_cvb.fh"
integer(kind=iwp), external :: mstackr_cvb
real(kind=wp), external :: tim0_cvb

cpu0 = tim0_cvb()
if (((.not. variat) .or. (nmcscf == 1)) .or. ((ip(3) >= 1) .and. ((.not. endvar) .or. (ip(6) >= 2)))) then
  cpu_prev = Zero
  n_applyt = 0
  n_applyh = 0
  n_hess = 0
  n_orbhess = 0
  n_cihess = 0
  n_2el = 0
  ibase0 = mstackr_cvb(0)
  call mfreer_cvb(ibase0)
  ibasemx = ibase0
else
  ibase0 = mstackr_cvb(0)
  call mfreer_cvb(ibase0)
  ibasemx = ibase0+memused
end if
n_iter = 0

return

end subroutine stat1_cvb
