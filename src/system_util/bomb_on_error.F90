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

function bomb_on_error() result(rc)

implicit none
logical :: rc
character(len=16) :: bomb, env

bomb = ' '
env = 'MOLCAS_BOMB'
call getenvf(env,bomb)
if (bomb(1:1) == 'Y' .or. bomb(1:1) == 'y' .or. bomb(1:1) == '1') then
  rc = .true.
else
  rc = .false.
end if

end function bomb_on_error
