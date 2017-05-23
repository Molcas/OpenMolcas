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
* Copyright (C) 2012, Thomas Bondo Pedersen                            *
************************************************************************
C#define _DEBUG_
      Subroutine ThisIsRestrictedCode(developer_name,code_name,Abort)
C
C     Thomas Bondo Pedersen, December 2012.
C
C     An extension of OnlyIMayUseIt, this routine allows you to
C     identify the code portion which is not for others to use.
C     Based on (and using) OnlyIMayUseIt by V. Veryazov.
C
C     If Abort: stop the execution
C
      Implicit None
      Character*(*) developer_name
      Character*(*) code_name
      Logical Abort
#include "warnings.fh"

      Character*256 val

      val=' '
      Call GetEnvf('MOLCAS_ISDEV',val)
      If (val.eq.'PRODUCTION') Return

#if defined (_DEBUG_)
      Call OnlyIMayUseIt(developer_name)
      Write(6,'(A,A)') '>>>>> Restricted code: ',code_name
      Write(6,'(A,A,//)') '>>>>> Contact ',developer_name
      If (Abort) Then
         If (val.eq.' ' .or. val.ne.developer_name) Then
            Call xQuit(_RC_GENERAL_ERROR_)
         End If
      End If
      Call xFlush(6)
#else
      If (val.eq.' ' .or. val.ne.developer_name) Then
         Call OnlyIMayUseIt(developer_name)
         Write(6,'(A,A,//)') '>>>>> Restricted code: ',code_name
         If (Abort) Call xQuit(_RC_GENERAL_ERROR_)
         Call xFlush(6)
      End If
#endif

      End
