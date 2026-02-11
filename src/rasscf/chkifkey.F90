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
      Subroutine ChkIfKey()
      use definitions, only: iwp
      use input_ras, only: LUInput, nKeys, CMD
      Implicit None
! ------------------------------------------------------------
! Check if the next item on luinput is a string that starts with
! a keyword. Print warning else.
! ------------------------------------------------------------
      Character(LEN=4) Command
      Character(LEN=180)  Line
      Integer(kind=iwp) iCmd
      Read(LUInput,*) Line
      Command=Line(1:4)
      Call UpCase(Command)
      Do iCmd=1,NKeys
        If ( Command.eq.Cmd(iCmd) ) Return
      End Do
      write(6,*)' ****************************************************'
      write(6,*)' ChkIfKey Warning: The following line seems intended'
      write(6,*)' to give some keyword input, but was not recognized:'
      write(6,*)' '''//line(1:32)//''''
      write(6,*)' Spelling or syntactic mistake? Ignored!'
      write(6,*)' ****************************************************'
      End Subroutine ChkIfKey
