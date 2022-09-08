!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1989-1992, Roland Lindh                                *
!               1990, IBM                                              *
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine Request_MCLR_Run(Run_MCLR,ireturn,iPrint)

use Definitions, only: iwp, u6

implicit none
logical(kind=iwp), intent(in) :: Run_MCLR
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp), intent(in) :: iPrint
#include "warnings.h"
integer(kind=iwp) :: LuInput
character(len=16) :: StdIn
integer(kind=iwp), external :: IsFreeUnit

if (Run_MCLR) then

  ! McKinley will automatically generate the input for MCLR
  ! and signal to AUTO (iRC=2) to run the input file Stdin.x.

  if (iPrint >= 6) then
    write(u6,*)
    write(u6,*) ' McKinley requests the MCLR module to be executed!'
    write(u6,*)
  end if

  LuInput = IsFreeUnit(11)
  call StdIn_Name(StdIn)
  call Molcas_Open(LuInput,StdIn)
  write(LuInput,'(A)') ' &MCLR &End'
  write(LuInput,'(A)') 'End of Input'
  close(LuInput)
  ireturn = _RC_INVOKED_OTHER_MODULE_
else
  ireturn = _RC_ALL_IS_WELL_
end if

return

end subroutine Request_MCLR_Run
