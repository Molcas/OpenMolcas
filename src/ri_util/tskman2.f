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
      Subroutine Init_Tsk2(id,mTask,jOpt,List)
#include "WrkSpc.fh"
      Integer List(*)  ! either nTask or 0 long
      Common /Tsk2/ iRsv,iOpt,nTask
*
      nTask=mTask
      iOpt=jOpt
      If (iOpt.eq.0) Then
         Call Init_Tsk(id,nTask)
      Else If (iOpt.eq.1) Then
         Call GetMem('TskList','Allo','Inte',id,nTask)
         Call ICopy(nTask,List,1,iWork(id),1)
         iRsv=0
      Else
         Call WarningMessage(2,'Error in Init_Tsk2')
         Write (6,*) 'Init_Tsk2: illegal iOpt value!'
         Call Abend()
      End If
*
      Return
      End
      Logical Function Rsv_Tsk2(id,kls)
      External Rsv_Tsk
#include "WrkSpc.fh"
      Common /Tsk2/ iRsv,iOpt,nTask
      Logical Rsv_Tsk
*
      If (iOpt.eq.0) Then
         Rsv_Tsk2=Rsv_Tsk(id,kls)
      Else If (iOpt.eq.1) Then
         Rsv_Tsk2=.True.
         If (iRsv+1.gt.nTask) Then
            Rsv_Tsk2=.False.
         Else
            kls=iWork(id+iRsv)
            iRsv=iRsv+1
            If (kls.le.0) Rsv_Tsk2=.False.
            If (kls.gt.nTask) Rsv_Tsk2=.False.
         End If
      Else
         Rsv_Tsk2=.False.
         Call WarningMessage(2,'Error in Rsv_Tsk2')
         Write (6,*) 'Rsv_Tsk2: illegal iOpt value!'
         Call Abend()
      End If
*
      Return
      End
      Subroutine Free_Tsk2(id)
      Common /Tsk2/ iRsv,iOpt,nTask
*
      If (iOpt.eq.0) Then
         Call Free_Tsk(id)
      Else If (iOpt.eq.1) Then
         Call GetMem('TskList','Free','Inte',id,nTask)
         nTask=0
      Else
         Call WarningMessage(2,'Error in Free_Tsk2')
         Write (6,*) 'Free_Tsk2: illegal iOpt value!'
         Call Abend()
      End If
      iOpt=-1
*
      Return
      End
