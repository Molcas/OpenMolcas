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

subroutine ChkIfKey()
! ------------------------------------------------------------
! Check if the next item on luinput is a string that starts with
! a keyword. Print warning else.
! ------------------------------------------------------------

use input_ras, only: LUInput, nKeys, CMD
use Definitions, only: iwp

implicit none
character(len=4) Command
character(len=180) Line
integer(kind=iwp) iCmd

read(LUInput,*) Line
Command = Line(1:4)
call UpCase(Command)
do iCmd=1,NKeys
  if (Command == Cmd(iCmd)) return
end do
write(6,*) ' ****************************************************'
write(6,*) ' ChkIfKey Warning: The following line seems intended'
write(6,*) ' to give some keyword input, but was not recognized:'
write(6,*) ' '''//line(1:32)//''''
write(6,*) ' Spelling or syntactic mistake? Ignored!'
write(6,*) ' ****************************************************'

end subroutine ChkIfKey
