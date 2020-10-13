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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Cho_TrcIdl_Update(IAmIdle)
C
C     Thomas Bondo Pedersen, May 2010.
C
C     Update array for tracing idle processors
C
      Implicit None
      Logical IAmIdle
#include "choptr2.fh"
#include "cho_para_info.fh"
#include "WrkSpc.fh"
#if defined (_DEBUGPRINT_)
#include "cholesky.fh"
#endif

#if defined (_DEBUGPRINT_)
      If (l_Idle.lt.1 .or. .not.Trace_Idle) Then
         Write(LuPri,'(A)')
     &   'Cho_TrcIdl_Update should not be called in this run!'
         Write(LuPri,*) 'Trace_Idle=',Trace_Idle
         Write(LuPri,'(A,2I10)')
     &   'ip_Idle,l_Idle=',ip_Idle,l_Idle
         Call Cho_Quit('Illegal call to Cho_TrcIdl_Update',103)
      End If
#endif

      If (IAmIdle) Then
         If (Cho_Real_Par) Then
            iWork(ip_Idle+myRank)=iWork(ip_Idle+myRank)+1
         Else
            iWork(ip_Idle)=iWork(ip_Idle)+1
         End If
      End If

      End
