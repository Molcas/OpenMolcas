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

integer IsStructure
character*256 value
character*100 Get_SuperName
external Get_SuperName

value = ' '
call getenvf('MOLCAS_STRUCTURE',value)
IsStructure = 0
if (value == '1') IsStructure = 1
if (Get_SuperName() == 'last_energy') IsStructure = 1

return

end function IsStructure
