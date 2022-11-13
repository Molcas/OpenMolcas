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

subroutine init_run_use()

use RunFile_data, only: i_run_CA_used, i_run_DA_used, i_run_DS_used, i_run_IA_used, i_run_IS_used

implicit none

i_run_CA_used(:) = 0
i_run_DA_used(:) = 0
i_run_DS_used(:) = 0
i_run_IA_used(:) = 0
i_run_IS_used(:) = 0

return

end subroutine init_run_use
