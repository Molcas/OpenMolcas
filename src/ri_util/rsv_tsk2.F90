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

logical function Rsv_Tsk2(id,kls)

use Tsk2

#include "stdalloc.fh"
logical, external :: Rsv_Tsk

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
  write(6,*) 'Rsv_Tsk2: illegal iOpt value!'
  call Abend()
end if

return

end function Rsv_Tsk2
