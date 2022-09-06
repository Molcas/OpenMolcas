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

subroutine Init_Tsk2(id,mTask,jOpt,List)

use Tsk2

#include "stdalloc.fh"
integer List(*)  ! either nTask or 0 long

nTask = mTask
iOpt = jOpt
if (iOpt == 0) then
  call Init_Tsk(id,nTask)
else if (iOpt == 1) then
  call mma_allocate(TskList,nTask,Label='TskList')
  TskList(1:nTask) = List(1:nTask)
  id = 0
  iRsv = 1
else
  call WarningMessage(2,'Error in Init_Tsk2')
  write(6,*) 'Init_Tsk2: illegal iOpt value!'
  call Abend()
end if

return

end subroutine Init_Tsk2
