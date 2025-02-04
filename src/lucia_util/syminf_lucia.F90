!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine SYMINF_LUCIA(IPRNT)
! Information about number of symmetries

use lucia_data, only: PNTGRP, NIRREP
use Definitions, only: u6

implicit none
integer IPRNT

if (PNTGRP == 1) then
  ! ===
  ! D2h
  ! ===
  call ZSYM1(NIRREP,IPRNT)
else
  write(u6,*) ' You are too early, sorry'
  write(u6,*) ' Illegal PNTGRP in SYMINF ',PNTGRP
  !stop 11
  call SYSABENDMSG('lucia_util/syminf','Internal error','')
end if

end subroutine SYMINF_LUCIA
