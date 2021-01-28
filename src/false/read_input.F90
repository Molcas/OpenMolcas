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
! Copyright (C) 2020, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Read_Input()

use False_Global, only: Run_Command
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: LU, EOF
character(len=4) :: Key
character(len=180) :: Line
integer(kind=iwp), external :: IsFreeUnit
character(len=180), external :: Get_Ln

LU = IsFreeUnit(11)
call SpoolInp(LU)
rewind(LU)
call RdNLst(LU,'false')

EOF = 0
do
  Line = Get_Ln(LU)
  if (EOF < 0) exit
  Key = Line(1:len(Key))
  call UpCase(Key)
  if (Key(1:3) == 'END') then
    exit
  else if (Key == 'RUN') then
    Line = Get_Ln(LU)
    if (EOF /= 0) then
      call WarningMessage(2,'Error reading RUN keyword.')
      call Quit_OnUserError()
    end if
    Run_Command = Line
  else
    call WarningMessage(2,'Error in FALSE input: Unknown keyword '//trim(key))
    call Quit_OnUserError()
  end if
end do

close(LU)

return

end subroutine Read_Input
