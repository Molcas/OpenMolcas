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

subroutine Cho_Free_Tsk_(ID)

use Para_Info, only: Is_Real_Par, Set_Do_Parallel
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ID
logical(kind=iwp) :: Parallel_Mode

Parallel_Mode = Is_Real_Par()

call Set_Do_Parallel(.false.)
call Free_Tsk(ID)
call Set_Do_Parallel(Parallel_Mode)

end subroutine Cho_Free_Tsk_
