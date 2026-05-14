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

function IsStructure()
!  This function return 0 if module is running outside structure
!                       1 if module is running inside structure
!                      -1 if module is running outside structure but there has been a structure loop earlier

use UnixInfo, only: SuperName
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: IsStructure
character(len=8) :: level, Val

Val = ' '
call GetEnvF('MOLCAS_STRUCTURE',Val)
IsStructure = 0
if (Val == '1') then
  IsStructure = -1
  call getenvf('EMIL_InLoop',level)
  if (level == '') level = '0'
  if (level(1:1) /= '0') IsStructure = 1
end if
if (SuperName == 'last_energy') IsStructure = 1

return

end function IsStructure
