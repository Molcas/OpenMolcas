!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

function Rsv_GTList(TskLw,TskHi,iOpt,NewBatch)

use definitions, only: iwp, wp
use Para_Info, only: nProcs, Is_Real_Par
use TList_Mod, only: iStrt_TList, iTCnSt, iTskCan, PQ
#ifdef _MOLCAS_MPP_
use TList_Mod, only: igaTsk, TskL, nTasks, iEnd_TList, mTasks, TskM
#endif
use Constants, only: One

implicit none
logical(kind=iwp) Rsv_GTList
logical(kind=iwp), intent(out) :: NewBatch
real(kind=wp), intent(out) :: TskLw, TskHi
integer(kind=iwp), intent(in) :: iOpt
#ifdef _MOLCAS_MPP_
integer(kind=iwp), external :: RsvTsk
integer(kind=iwp) MyTask
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
Rsv_GTList = .false.
if ((.not. Is_Real_Par()) .or. (nProcs == 1)) then
  if (iTCnSt == 1) then
    iTCnSt = iTCnSt+1
    iTskCan = iTskCan+1
    Rsv_GTList = .true.
    TskLw = One
    TskHi = PQ
    iStrt_TList = 1
    if (iOpt == 0) then
      NewBatch = .true.
    else
      NewBatch = .false.
    end if
  end if
#ifdef _MOLCAS_MPP_
else
  !
  !---- NewBatch is true if
  !      1) Batch never processed by the node before
  !      2) The sequence of executed batches is broken

  if (iOpt == 0) then
    MyTask = RsvTsk(igaTsk,TskL,nTasks,nTasks,iTCnST,iStrt_TList,iEnd_TList)
    NewBatch = .true.
  else if (iOpt == 1) then
    MyTask = RsvTsk(igaTsk,TskL,nTasks,mTasks,iTCnST,iStrt_TList,iEnd_TList)
    NewBatch = iStrt_TList > mTasks
  else if (iOpt == 2) then
    MyTask = RsvTsk(igaTsk,TskL,nTasks,nTasks,iTCnST,iStrt_TList,iEnd_TList)
    NewBatch = (iStrt_TList > mTasks) .or. (iStrt_TList /= iTCnST)
  else
    MyTask = 0
    write(6,*) 'Rsv_GTList: Invalid option:',iOpt
    call Abend()
  end if
  if (MyTask >= 1) then
    Rsv_GTList = .true.
    TskLw = TskM(1,MyTask)
    TskHi = TskM(2,MyTask)
    iTCnSt = iTCnSt+1
    iTskCan = iTskCan+1
  end if
#endif
end if
!                                                                      *
!***********************************************************************
!                                                                      *
end function Rsv_GTList
