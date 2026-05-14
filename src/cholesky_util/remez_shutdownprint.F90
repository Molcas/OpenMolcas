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

subroutine Remez_ShutDownPrint(print_to_molcas_log)

use ReMez_mod, only: IW
use Definitions, only: iwp

implicit none
logical(kind=iwp), intent(in) :: print_to_molcas_log

if ((.not. print_to_molcas_log) .and. (IW > 0)) then
  close(IW)
  IW = -999
end if

end subroutine Remez_ShutDownPrint
