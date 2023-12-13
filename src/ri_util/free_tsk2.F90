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

subroutine Free_Tsk2(id)

use RI_glob, only: iOpt, nTask, TskList
use stdalloc, only: mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: id

if (iOpt == 0) then
  call Free_Tsk(id)
else if (iOpt == 1) then
  call mma_deallocate(TskList)
  nTask = 0
else
  call WarningMessage(2,'Error in Free_Tsk2')
  write(u6,*) 'Free_Tsk2: illegal iOpt value!'
  call Abend()
end if
iOpt = -1

return

end subroutine Free_Tsk2
