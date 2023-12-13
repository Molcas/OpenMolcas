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
! Copyright (C) 2007, Ten-no Research Group                            *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Remez_SetupPrint(print_to_molcas_log)

use ReMez_mod, only: IW
use Definitions, only: iwp, u6

implicit none
logical(kind=iwp), intent(in) :: print_to_molcas_log
integer(kind=iwp) :: ini
integer(kind=iwp), external :: isFreeUnit

if (print_to_molcas_log) then
  IW = u6
else
  ini = 7
  IW = isFreeUnit(ini)
  call Molcas_Open(IW,'REMEZ')
end if

end subroutine Remez_SetupPrint
