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
      use Tsk2
#include "stdalloc.fh"
      Integer List(*)  ! either nTask or 0 long
*
      nTask=mTask
      iOpt=jOpt
      If (iOpt.eq.0) Then
         Call Init_Tsk(id,nTask)
      Else If (iOpt.eq.1) Then
         Call mma_allocate(TskList,nTask,Label='TskList')
         TskList(1:nTask) = List(1:nTask)
         id = 0
         iRsv=1
      Else
         Call WarningMessage(2,'Error in Init_Tsk2')
         Write (6,*) 'Init_Tsk2: illegal iOpt value!'
         Call Abend()
      End If
*
      Return
      End
      Logical Function Rsv_Tsk2(id,kls)
      use Tsk2
#include "stdalloc.fh"
      Logical, External :: Rsv_Tsk
*
      If (iOpt.eq.0) Then
         Rsv_Tsk2=Rsv_Tsk(id,kls)
      Else If (iOpt.eq.1) Then
         Rsv_Tsk2=.True.
         If (iRsv.gt.nTask) Then
            Rsv_Tsk2=.False.
         Else
            kls=TskList(iRsv)
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
      use Tsk2
#include "stdalloc.fh"
*
      If (iOpt.eq.0) Then
         Call Free_Tsk(id)
      Else If (iOpt.eq.1) Then
         Call mma_deallocate(TskList)
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
