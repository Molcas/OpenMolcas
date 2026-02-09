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
      use definitions, only: iwp, u6
      use TList_Mod, only: iStrt_TList,iEnd_TList,iTskCan,mTasks,
     &                     Not_Used,nTasks,PP_Status,QLast, TskL
      Use Para_Info, only: MyRank, nProcs, Is_Real_Par
      Implicit None
      Logical(kind=iwp):: Debug=.False.
      Integer(kind=iwp), Pointer:: TskList(:,:)
      Integer(kind=iwp) i, iE, iTsk
*
      If (Debug) Then
         If (PP_Status) Then
            Write (u6,*) 'Init_PPList: Active'
         Else
            Write (u6,*) 'Init_PPList: InActive'
         End If
      End If

      If (PP_Status) Return
      PP_Status=.True.
*
      iTskCan=0
      mTasks=0
      iStrt_TList=0
      iEnd_TList=nTasks+1
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return

      TskList(1:nTasks,1:2) => TskL(1:2*nTasks)
      TskList(:,1)=0
      Do iTsk = 0, nTasks-1
        TskList(1+iTsk,1)=MOD(iTsk+MyRank,nTasks)+1
      End Do
*
*---- Copy list in inverse order to ensure the proper order if
*     reinitiated at once.
*
      iE = nTasks-1
      TskList(:,2)=0
      Do i = 0, nTasks-1
         TskList(1+iE,2) = TskList(i+1,1)
         iE = iE - 1
      End Do
      nullify(TskList)
*
      QLast(1)=Not_Used
      QLast(2)=Not_Used
*
      End Subroutine Init_PPList

      Subroutine ReInit_PPList(Semi_Direct)
      use definitions, only: iwp, u6
      Use Para_Info, only: MyRank, nProcs
      use TList_Mod, only: TskL,iStrt_TList,iEnd_TList,iTskCan,mTasks,
     &                     Not_Used,nTasks,PP_Status,QLast
      Implicit None
      Logical(kind=iwp), intent(in):: Semi_Direct
      Logical(kind=iwp):: Debug=.False.
      Integer(kind=iwp), Pointer:: TskList(:,:)
      Integer(kind=iwp) i, iCount, iE
*
      If (Debug) Then
         If (PP_Status) Then
            Write (u6,*) 'ReInit_PPList: Active'
         Else
            Write (u6,*) 'ReInit_PPList: InActive'
         End If
      End If
      If (.NOT.PP_Status) Then
         Write (u6,*) 'ReInit_PPList: List is not active!'
         Call Abend()
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
         Call ICopy(mTasks,TskList(:,2),1,TskList(:,1),1)
*
*---- Now copy task indices of tasks which were not executed by this node.
*     Change the order so that the first task it the largest in the list.
*
         iE = mTasks
         iCount = 1
         Do i = mTasks, nTasks-1
            If (iCount .gt. myRank) Then
               TskList(1+i,1) = TskList(i+1,2)
            Else
               TskList(1+i,1) = TskList(iE+1,2)
               iE = iE - 1
               iCount = iCount + 1
            Endif
         Enddo
*
         nullify(TskList)
      End If
*
      iStrt_TList=0
      iEnd_TList=nTasks+1
*
      QLast(1)=Not_Used
      QLast(2)=Not_Used
*
      End Subroutine ReInit_PPList
*
      Subroutine Free_PPList()
      use TList_Mod, only: TskL, PP_Status
      Use Para_Info, only: nProcs, Is_Real_Par
      Use stdalloc, Only: mma_deallocate
      implicit none
*
      If (.NOT.Allocated(TskL)) Return
      PP_Status=.False.
*
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return
      Call mma_deallocate(TskL)
*
      End Subroutine Free_PPList
