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
      Subroutine Init_GTList()
      use TList_Mod, only: GT_Status,iTCnSt
#ifdef _MOLCAS_MPP_
      use TList_Mod, only: nTasks,igaTsk
#endif
      Use Para_Info, Only: nProcs, Is_Real_Par
      implicit none
*
      If (GT_Status) Return
      GT_Status=.True.
*
      iTCnSt = 1
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return
#ifdef _MOLCAS_MPP_
*     create global tasklist...
      Call GATskL(.TRUE.,nTasks,igaTsk)
#endif
*
      End Subroutine Init_GTList

      Subroutine ReInit_GTList()
      use TList_Mod, only: GT_Status,iTCnSt
#ifdef _MOLCAS_MPP_
      use TList_Mod, only: iGATsk
#endif

      Use Para_Info, Only: nProcs, Is_Real_Par
      implicit none
*
      If (.Not.GT_Status) Then
         Write (6,*) 'ReInit_GTList: List not active!'
         Call Abend()
      End If
      iTCnSt = 1
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return
#ifdef _MOLCAS_MPP_
*     initialize global tasklist...
      Call GATskL_Zero(igaTsk)
#endif
*
      End Subroutine ReInit_GTList

************************************************************************
* Rsv_GTList
************************************************************************
      Function Rsv_GTList(TskLw,TskHi,iOpt,NewBatch)
      use definitions, only: iwp, wp
      Use Para_Info, Only: nProcs, Is_Real_Par
      use TList_Mod, only: iStrt_TList,iTCnSt,iTskCan,PQ
#ifdef _MOLCAS_MPP_
      use TList_Mod, only: igaTsk,TskL,nTasks,iEnd_TList,mTasks,TskM
#endif
      use Constants, only: One
      Implicit None
      Logical(kind=iwp) Rsv_GTList
      Logical(kind=iwp), intent(out):: NewBatch
      real(kind=wp), intent(out):: TskLw,TskHi
      integer(kind=iwp), intent(in):: iOpt
#ifdef _MOLCAS_MPP_
      Integer(kind=iwp), External :: RsvTsk
      Integer(kind=iwp) MyTask
#endif
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
            MyTask=RsvTsk(igaTsk,TskL,nTasks,nTasks,iTCnST,
     &                    iStrt_TList,iEnd_TList)
            NewBatch = .True.
         Else If (iOpt.eq.1) Then
            MyTask=RsvTsk(igaTsk,TskL,nTasks,mTasks,iTCnST,
     &                    iStrt_TList,iEnd_TList)
            NewBatch = iStrt_TList.gt.mTasks
         Else If (iOpt.eq.2) Then
            MyTask=RsvTsk(igaTsk,TskL,nTasks,nTasks,iTCnST,
     &                    iStrt_TList,iEnd_TList)
            NewBatch = iStrt_TList.gt.mTasks .or.
     &                 iStrt_TList.ne.iTCnST
         Else
            MyTask=0
            Write (6,*) 'Rsv_GTList: Invalid option:',iOpt
            Call Abend()
         End If
         If (MyTask.ge.1) Then
            Rsv_GTList=.True.
            TskLw=TskM(1,MyTask)
            TskHi=TskM(2,MyTask)
            iTCnSt=iTCnSt+1
            iTskCan=iTskCan+1
         End If
#endif
      End If
*                                                                      *
************************************************************************
*                                                                      *
      End Function Rsv_GTList

************************************************************************
* Free_GTList
************************************************************************
      Subroutine Free_GTList()
      Use Para_Info, Only: nProcs, Is_Real_Par
      use TList_Mod, only: GT_Status,iTCnSt
#ifdef _MOLCAS_MPP_
      use TList_Mod, only: nTasks,igaTsk
#endif
      Implicit None
*
      If (.NOT.GT_Status) Return
      GT_Status=.False.
*
      iTCnSt = 1
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return
#ifdef _MOLCAS_MPP_
*     create global tasklist...
      Call GATskL(.False.,nTasks,igaTsk)
#endif
*
      End Subroutine Free_GTList
