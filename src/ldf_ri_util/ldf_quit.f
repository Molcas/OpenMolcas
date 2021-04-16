************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine LDF_Quit(irc)
C
C     Quit with special error codes:
C       irc=-1 => _RC_INPUT_ERROR_
C       irc=1  => _RC_INTERNAL_ERROR
C       other  => _RC_GENERAL_ERROR
C
      Implicit None
      Integer irc
#include "warnings.fh"

      If (irc.eq.-1) Then
         Call xQuit(_RC_INPUT_ERROR_)
      Else If (irc.eq.1) Then
         Call xQuit(_RC_INTERNAL_ERROR_)
      Else
         Call xQuit(_RC_GENERAL_ERROR_)
      End If

      End
