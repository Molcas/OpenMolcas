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

subroutine LDF_CheckConfig()

#ifdef _MOLCAS_MPP_
use Para_Info, only: nProcs, Is_Real_Par
#endif
use Definitions, only: u6

implicit none
#include "localdf.fh"

! Debug write of unconstrained coefficients:
! 1) makes no sense for unconstrained LDF => reset
! 2) not implemented in parallel => error
WriteUnconstrainedC = WriteUnconstrainedC .and. (LDF_Constraint /= -1)
#ifdef _MOLCAS_MPP_
if (WriteUnconstrainedC) then
  if ((nProcs > 1) .and. Is_Real_Par()) then
    call WarningMessage(2,'Write unconstrained coefficients not implemented in parallel')
    call Quit_OnUserError()
  end if
end if
#endif
! Using unique atom pairs is buggy, warn!
if (UseUniqueAtomPairs) then
  call WarningMessage(1,'WARNING: using unique atom pairs may cause erroneous results')
  call xFlush(u6)
end if

return

end subroutine LDF_CheckConfig
