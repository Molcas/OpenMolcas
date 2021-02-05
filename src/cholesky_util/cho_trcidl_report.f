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
      SubRoutine Cho_TrcIdl_Report()
C
C     Thomas Bondo Pedersen, May 2010.
C
C     Report idle status for all processors
C
      Use Para_Info, Only: nProcs
      use ChoArr, only: Idle
      Implicit None
#include "cholesky.fh"
#include "cho_para_info.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Integer i
      Integer nIdle
      Integer, Allocatable:: TIloc(:)

#if defined (_DEBUGPRINT_)
      If (.NOT.Allocated(Idle) .or. .not.Trace_Idle) Then
         Write(LuPri,'(A)')
     &   'Cho_TrcIdl_Report should not be called in this run!'
         Write(LuPri,*) 'Trace_Idle=',Trace_Idle
         Call Cho_Quit('Illegal call to Cho_TrcIdl_Report',103)
      End If
#endif

      If (Cho_Real_Par) Then
#if defined (_DEBUGPRINT_)
         l_Idle=0
         If (Allocated(Idle)) l_Idle = SIZE(Idle)
         If (l_Idle.lt.nProcs) Then
            Write(LuPri,'(A)')
     &      'Error detected in Cho_TrcIdl_Report: l_Idle < nProcs'
            Write(LuPri,*) 'Trace_Idle=',Trace_Idle
            Call Cho_Quit(
     &               'Cho_TrcIdle_Report: l_Idle not properly set!',103)
         End If
#endif
         Call mma_allocate(TIloc,[0,nProcs-1],Label='TILoc')
         Call iCopy(nProcs,Idle,1,TIloc,1)
         Call Cho_GAIGOp(TIloc,nProcs,'+')
         nIdle=0
         Do i=0,nProcs-1
            nIdle=nIdle+min(TIloc(i),1)
         End Do
         If (nIdle.eq.0) Then
            Write(LuPri,'(A)')
     &      'No idle procs to report'
         Else
            Write(LuPri,'(I4,A,I4,A,F7.2,A)')
     &      nIdle,' of',nProcs,' procs have been idle (',
     &      1.0d2*dble(nIdle)/dble(nProcs),' %)'
            Write(LuPri,'(A)')
     &      'List of idle procs:'
            Do i=0,nProcs-1
               If (TIloc(i).gt.0) Then
                  Write(LuPri,'(I4,A,I8,A)')
     &            i,' (Idle counter:',TIloc(i),')'
               End If
            End Do
         End If
         Call mma_deallocate(TIloc)
      Else
         If (Idle(1).eq.0) Then
            Write(LuPri,'(A)')
     &      'No idle procs to report!'
         Else
            Write(LuPri,'(A,I8,A)')
     &      'Proc 0 has been idle',Idle(1),' times'
         End If
      End If
      Call Cho_Flush(LuPri)

      End
