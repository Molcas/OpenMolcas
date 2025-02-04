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
! Copyright (C) 1991, Jeppe Olsen                                      *
!***********************************************************************

subroutine SYMCOM(ITASK,I1,I2,I12)
! Symmetries I1,I2,I12 are related as
! I1*I2 = I12
! IF (ITASK == 1) I2 and I12 are known, find I1
! IF (ITASK == 2) I1 and I12 are known, find I1
! IF (ITASK == 3) I1 and I2 are known, find I12
!
! If obtained symmetry I1 or I2 is outside bounds,
! zero is returned.
!
! Jeppe Olsen, Spring of 1991
!
! ==============
! Driver routine
! ==============

use lucia_data, only: PNTGRP
use Definitions, only: u6

implicit none
integer ITASK, I1, I2, I12

if (PNTGRP == 1) then
  call SYMCM1(ITASK,I1,I2,I12)
else
  write(u6,*) ' PNTGRP parameter out of bounds ',PNTGRP
  write(u6,*) ' Enforced stop in SYMCOM'
  call SYSABENDMSG('lucia_util/symcom','Internal error','')
end if

end subroutine SYMCOM
