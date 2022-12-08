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

function Rsv_Tsk2(id,kls)

use RI_glob, only: iOpt, iRsv, nTask, TskList
use Definitions, only: iwp, u6

implicit none
logical(kind=iwp) :: Rsv_Tsk2
integer(kind=iwp), intent(in) :: id
integer(kind=iwp), intent(out) :: kls
logical(kind=iwp), external :: Rsv_Tsk

if (iOpt == 0) then
  Rsv_Tsk2 = Rsv_Tsk(id,kls)
else if (iOpt == 1) then
  Rsv_Tsk2 = .true.
  if (iRsv > nTask) then
    Rsv_Tsk2 = .false.
  else
    kls = TskList(iRsv)
    iRsv = iRsv+1
    if (kls <= 0) Rsv_Tsk2 = .false.
    if (kls > nTask) Rsv_Tsk2 = .false.
  end if
else
  Rsv_Tsk2 = .false.
  call WarningMessage(2,'Error in Rsv_Tsk2')
  write(u6,*) 'Rsv_Tsk2: illegal iOpt value!'
  call Abend()
end if

return

end function Rsv_Tsk2
