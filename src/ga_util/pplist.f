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
      Subroutine Init_PPList()
      use TList_Mod
      Use Para_Info, only: MyRank, nProcs, Is_Real_Par
      Implicit Real*8 (a-h,o-z)
#include "status.fh"
      Logical:: Debug=.False.
      Integer, Pointer:: TskList(:,:)
*
      If (Debug) Then
         If (PP_Status.eq.Active) Then
            Write (6,*) 'Init_PPList: Active'
         Else If (PP_Status.eq.InActive) Then
            Write (6,*) 'Init_PPList: Active'
         Else
            Write (6,*) 'Init_PPList: undefined!'
         End If
      End If
      If (PP_Status.eq.Active) Return
      PP_Status=Active
*
      Write (6,*) 'Init_PPList'
      Call xflush(6)
      iTskCan=0
      mTasks=0
      iStrt_TList=0
      iEnd_TList=nTasks+1
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return

      TskList(1:nTasks,1:2) => TskL(1:2*nTasks)
*     call izero(TskL(1:nTasks),nTasks)
      TskList(:,1)=0
      Do iTsk = 0, nTasks-1
*       TskL(1+iTsk)=MOD(iTsk+MyRank,nTasks)+1
        TskList(1+iTsk,1)=MOD(iTsk+MyRank,nTasks)+1
      End Do
c     Write (*,*) (TskL(iTsk),iTsk = 1, nTasks)
*
*---- Copy list in inverse order to ensure the proper order if
*     reinitiated at once.
*
      iE = nTasks-1
*     call izero(TskL(1+nTasks:2*nTasks),nTasks)
      TskList(:,2)=0
      Do i = 0, nTasks-1
*        TskL(1+nTasks+iE) = TskL(i+1)
         TskList(1+iE,2) = TskList(i+1,1)
         iE = iE - 1
      End Do
      TskList => Null()
*
      QLast(1)=Not_Used
      QLast(2)=Not_Used
*
      Return
      End
      Subroutine ReInit_PPList(Semi_Direct)
      Use Para_Info, only: MyRank, nProcs
      use TList_Mod
      Implicit Real*8 (a-h,o-z)
#include "status.fh"
      Logical Semi_Direct
      Logical:: Debug=.False.
      Integer, Pointer:: TskList(:,:)
*
      Write (6,*) 'ReInit_PPList'
      Call xflush(6)
*
      If (Debug) Then
         If (PP_Status.eq.Active) Then
            Write (6,*) 'ReInit_PPList: Active'
         Else If (PP_Status.eq.InActive) Then
            Write (6,*) 'ReInit_PPList: Active'
         Else
            Write (6,*) 'ReInit_PPList: undefined!'
         End If
      End If
      If (PP_Status.ne.Active) Then
         Write (6,*) 'ReInit_PPList: List is not active!'
         Call Abend
      End If
      iTskCan=0
      mTasks=iStrt_TList
      If (nProcs.eq.1) Then
         iStrt_TList=0
         iEnd_TList=nTasks+1
         Return
      End If
*
      If (Semi_Direct) Then
         TskList(1:nTasks,1:2) => TskL(1:2*nTasks)
*
*---- Copy first the task indices of tasks that were exectuted
*        Call ICopy(mTasks,TskL(1+nTasks:2*nTasks),1,TskL(1:mTasks),1)
         Call ICopy(mTasks,TskList(:,2),1,TskList(:,1),1)
c     Write (*,*) 'mTasks=',mTasks
*
*---- Now copy task indices of tasks which were not executed by this node.
*     Change the order so that the first task it the largest in the list.
*
         iE = mTasks
         iCount = 1
         Do i = mTasks, nTasks-1
            If (iCount .gt. myRank) Then
*              TskL(1+i) = TskL(i+1+nTasks)
               TskList(1+i,1) = TskList(i+1,2)
            Else
*              TskL(1+i) = TskL(iE+1+nTasks)
               TskList(1+i,1) = TskList(iE+1,2)
               iE = iE - 1
               iCount = iCount + 1
            Endif
         Enddo

c        Write (*,*) (TskL(iTsk),iTsk = 1, nTasks)
c        Write (*,*) 'mTasks=',mTasks
*
         TskList => Null()
      End If
*
      iStrt_TList=0
      iEnd_TList=nTasks+1
*
      QLast(1)=Not_Used
      QLast(2)=Not_Used
*
      Return
      End
*
      Subroutine Free_PPList()
      use TList_Mod
      Use Para_Info, only: nProcs, Is_Real_Par
#include "status.fh"
#include "stdalloc.fh"
*
*     If (PP_Status.ne.Active) Return
      If (.NOT.Allocated(TskL)) Return
      PP_Status=Inactive
*
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return
      Call mma_deallocate(TskL)
*
      Return
      End
