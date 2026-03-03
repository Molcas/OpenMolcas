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

subroutine Init_PPList()

use definitions, only: iwp, u6
use TList_Mod, only: iStrt_TList, iEnd_TList, iTskCan, mTasks, Not_Used, nTasks, PP_Status, QLast, TskL
use Para_Info, only: MyRank, nProcs, Is_Real_Par

implicit none
logical(kind=iwp) :: Debug = .false.
integer(kind=iwp), pointer :: TskList(:,:)
integer(kind=iwp) i, iE, iTsk

if (Debug) then
  if (PP_Status) then
    write(u6,*) 'Init_PPList: Active'
  else
    write(u6,*) 'Init_PPList: InActive'
  end if
end if

if (PP_Status) return
PP_Status = .true.

iTskCan = 0
mTasks = 0
iStrt_TList = 0
iEnd_TList = nTasks+1
if ((.not. Is_Real_Par()) .or. (nProcs == 1)) return

TskList(1:nTasks,1:2) => TskL(1:2*nTasks)
TskList(:,1) = 0
do iTsk=0,nTasks-1
  TskList(1+iTsk,1) = mod(iTsk+MyRank,nTasks)+1
end do

! Copy list in inverse order to ensure the proper order if
! reinitiated at once.

iE = nTasks-1
TskList(:,2) = 0
do i=0,nTasks-1
  TskList(1+iE,2) = TskList(i+1,1)
  iE = iE-1
end do
nullify(TskList)

QLast(1) = Not_Used
QLast(2) = Not_Used

end subroutine Init_PPList
