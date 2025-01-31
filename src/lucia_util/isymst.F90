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

integer function ISYMST(STRING,NEL)
! Master routine for symmetry of string

use lucia_data, only: PNTGRP

implicit none
integer NEL
integer STRING(*)
integer, external :: ISYMS1

if (PNTGRP == 1) then
  ! D2h
  ISYMST = ISYMS1(STRING,NEL)
else
  write(6,*) ' Sorry PNTGRP option not programmed ',PNTGRP
  write(6,*) ' Enforced stop in ISYMST'
  !stop 5
  call SYSABENDMSG('lucia_util/isymst','Internal error','')
  ISYMST = -9999
end if

end function ISYMST
