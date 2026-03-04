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

subroutine ReInit_PPList(Semi_Direct)

use Para_Info, only: MyRank, nProcs
use TList_Mod, only: iEnd_TList, iStrt_TList, iTskCan, mTasks, Not_Used, nTasks, PP_Status, QLast, TskL
use Definitions, only: iwp, u6

implicit none
logical(kind=iwp), intent(in) :: Semi_Direct
integer(kind=iwp) i, iCount, iE
integer(kind=iwp), pointer :: TskList(:,:)
logical(kind=iwp), parameter :: Debug = .false.

if (Debug) then
  if (PP_Status) then
    write(u6,*) 'ReInit_PPList: Active'
  else
    write(u6,*) 'ReInit_PPList: InActive'
  end if
end if
if (.not. PP_Status) then
  write(u6,*) 'ReInit_PPList: List is not active!'
  call Abend()
end if
iTskCan = 0
mTasks = iStrt_TList
if (nProcs == 1) then
  iStrt_TList = 0
  iEnd_TList = nTasks+1
else
  if (Semi_Direct) then
    TskList(1:nTasks,1:2) => TskL(1:2*nTasks)

    ! Copy first the task indices of tasks that were exectuted
    call ICopy(mTasks,TskList(:,2),1,TskList(:,1),1)

    ! Now copy task indices of tasks which were not executed by this node.
    ! Change the order so that the first task it the largest in the list.

    iE = mTasks
    iCount = 1
    do i=mTasks,nTasks-1
      if (iCount > myRank) then
        TskList(1+i,1) = TskList(i+1,2)
      else
        TskList(1+i,1) = TskList(iE+1,2)
        iE = iE-1
        iCount = iCount+1
      end if
    end do

    nullify(TskList)
  end if

  iStrt_TList = 0
  iEnd_TList = nTasks+1

  QLast(1) = Not_Used
  QLast(2) = Not_Used
end if

end subroutine ReInit_PPList
