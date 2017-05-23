************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2013, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine RPA_Warn(Level,Message)
C
C     Thomas Bondo Pedersen (CTCC,UiO), July 2013.
C
C     Level <= 1 --- issue warning with message Message and return.
C     Level >= 2 --- issue warning with message Message and quit.
C                    Error codes used in xQuit:
C                    Level=2: _RC_INPUT_ERROR_
C                    Level=3: _RC_INTERNAL_ERROR_
C                    Level>3: _RC_GENERAL_ERROR_
C
      Implicit None
      Integer Level
      Character*(*) Message
#include "warnings.fh"

      Integer iLevel
      Integer rc

      If (Level.le.1) Then
         iLevel=max(Level,0)
         rc=0
      Else
         iLevel=2
         If (Level.eq.2) Then
            rc=_RC_INPUT_ERROR_
         Else If (Level.eq.3) Then
            rc=_RC_INTERNAL_ERROR_
         Else
            rc=_RC_GENERAL_ERROR_
         End If
      End If
      Call WarningMessage(iLevel,Message)
      Call xFlush(6)
      If (rc.ne.0) Then
         Call xQuit(rc)
      End If

      End
