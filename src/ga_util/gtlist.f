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
************************************************************************
* Init_GTList
************************************************************************
      Subroutine Init_GTList
      use TList_Mod
      Use Para_Info, Only: nProcs, Is_Real_Par
#include "status.fh"
*
      If (GT_Status.eq.Active) Return
      GT_Status=Active
*
      iTCnSt = 1
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return
#ifdef _MOLCAS_MPP_
*     create global tasklist...
      Call GATskL(.TRUE.,nTasks,igaTsk)
#endif
*
      Return
      End
      Subroutine ReInit_GTList
      use TList_Mod
      Use Para_Info, Only: nProcs, Is_Real_Par
#include "status.fh"
*
      If (GT_Status.ne.Active) Then
         Write (6,*) 'ReInit_GTList: List not active!'
         Call Abend
      End If
      iTCnSt = 1
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return
#ifdef _MOLCAS_MPP_
*     initialize global tasklist...
      Call GATskL_Zero(igaTsk)
#endif
*
      Return
      End

************************************************************************
* Rsv_GTList
************************************************************************
      Logical Function Rsv_GTList(TskLw,TskHi,iOpt,NewBatch)
      Use Para_Info, Only: nProcs, Is_Real_Par
      use TList_Mod
      Implicit Real*8 (a-h,o-z)
#ifdef _MOLCAS_MPP_
      External RsvTsk
      Integer RsvTsk
#endif
      Logical NewBatch
#include "tlist.fh"
#include "WrkSpc.fh"
#include "real.fh"
*                                                                      *
************************************************************************
*                                                                      *
      Rsv_GTList=.False.
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Then
         If (iTCnSt.eq.1) Then
            iTCnSt=iTCnSt+1
            iTskCan=iTskCan+1
            Rsv_GTList=.True.
            TskLw=One
            TskHi=PQ
            iStrt_TList=1
            If (iOpt.eq.0) Then
               NewBatch = .True.
            Else
               NewBatch = .False.
            End If
         End If
#ifdef _MOLCAS_MPP_
      Else
*
*------ NewBatch is true if
*        1) Batch never processed by the node before
*        2) The sequence of executed batches is broken
*
         If (iOpt.eq.0) Then
            MyTask=RsvTsk(igaTsk,iWork(ipTskL),nTasks,nTasks,iTCnST,
     &                    iStrt_TList,iEnd_TList)
            NewBatch = .True.
         Else If (iOpt.eq.1) Then
            MyTask=RsvTsk(igaTsk,iWork(ipTskL),nTasks,mTasks,iTCnST,
     &                    iStrt_TList,iEnd_TList)
            NewBatch = iStrt_TList.gt.mTasks
         Else If (iOpt.eq.2) Then
            MyTask=RsvTsk(igaTsk,iWork(ipTskL),nTasks,nTasks,iTCnST,
     &                    iStrt_TList,iEnd_TList)
            NewBatch = iStrt_TList.gt.mTasks .or.
     &                 iStrt_TList.ne.iTCnST
         Else
            MyTask=0
            Write (6,*) 'Rsv_GTList: Invalid option:',iOpt
            Call Abend
         End If
         If (MyTask.ge.1) Then
            Rsv_GTList=.True.
            TskLw=Work(ipTskM+2*(MyTask-1))
            TskHi=Work(ipTskM+2*(MyTask-1)+1)
            iTCnSt=iTCnSt+1
            iTskCan=iTskCan+1
         End If
#endif
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End

************************************************************************
* Free_GTList
************************************************************************
      Subroutine Free_GTList
      Use Para_Info, Only: nProcs, Is_Real_Par
      use TList_Mod
#include "status.fh"
*
      If (GT_Status.ne.Active) Return
      GT_Status=Inactive
*
      iTCnSt = 1
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return
#ifdef _MOLCAS_MPP_
*     create global tasklist...
      Call GATskL(.False.,nTasks,igaTsk)
#endif
*
      Return
      End
