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

function ISYMST(STRING,NEL)
! Master routine for symmetry of string

use lucia_data, only: PNTGRP
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: ISYMST
integer(kind=iwp) :: STRING(*), NEL
integer(kind=iwp), external :: ISYMS1

if (PNTGRP == 1) then
  ! D2h
  ISYMST = ISYMS1(STRING,NEL)
else
  write(u6,*) ' Sorry PNTGRP option not programmed ',PNTGRP
  write(u6,*) ' Enforced stop in ISYMST'
  !stop 5
  call SYSABENDMSG('lucia_util/isymst','Internal error','')
  ISYMST = -9999
end if

end function ISYMST
