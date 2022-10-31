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

subroutine IniStat()

use Para_Info, only: nProcs
use Definitions, only: u6

implicit none
#include "WrkSpc.fh"
#include "timtra.fh"

if (nfld_stat == 0) return
if (nfld_stat > nfldmax) then
  call WarningMessage(2,'Too many fields in IniStat')
  write(u6,*) 'nfld_stat:',nfld_stat
  call Abend()
end if
call GetMem('iGAStat','Allo','Real',iGAStat,nProcs*nfld_stat)
call Fzero(Work(iGAStat),nProcs*nfld_stat)

return

end subroutine IniStat
