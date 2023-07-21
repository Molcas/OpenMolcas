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

logical function Cho_Rsv_Tsk_(ID,i)

use Para_Info, only: Is_Real_Par, Set_Do_Parallel

implicit none
integer ID, i
#include "cho_para_info.fh"
logical Parallel_Mode
logical Rsv_Tsk
external Rsv_Tsk

Parallel_Mode = Is_Real_Par()

call Set_Do_Parallel(.false.)
Cho_Rsv_Tsk_ = Rsv_Tsk(ID,i)
call Set_Do_Parallel(Parallel_Mode)

end function Cho_Rsv_Tsk_
