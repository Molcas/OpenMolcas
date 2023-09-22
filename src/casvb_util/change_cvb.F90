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

subroutine change_cvb()

use Constants, only: Ten
use Definitions, only: iwp

implicit none
#include "main_cvb.fh"
#include "files_cvb.fh"
integer(kind=iwp) :: kbasisp
logical(kind=iwp), external :: chpcmp_cvb, lchpcmp_cvb, & ! ... Change of dimensioning variables ...
                               up2date_cvb ! ... Make: up to date? ...

! General settings:
proj = projsym

! Determine changes to memory assignment:
call chpcmp0_cvb()

call change1_cvb()
call change2_cvb()
call change3_cvb()
call change4_cvb()
call change5_cvb()
call change6_cvb()
call change7_cvb()

! Finished determining reassignment of memory, now determine
! what information must be evaluated/read again:

call chpcmp2_cvb(kbasis,kbasisp)
if (up2date_cvb('GUESS') .and. (kbasiscvb /= kbasis)) call touch_cvb('TRNSPN')
call symchk_cvb()
if (chpcmp_cvb(nint(strtint*Ten))) call touch_cvb('RDINT')

! Redo CIVB if definition has changed (from CVB or CIVECP):
if (lchpcmp_cvb(projcas)) then
  call setcnt2_cvb(2,0)
  call setcnt2_cvb(3,0)
  call setcnt2_cvb(4,0)
  call setcnt2_cvb(5,0)
end if

return

end subroutine change_cvb
