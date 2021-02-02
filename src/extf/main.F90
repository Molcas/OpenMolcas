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

program Main

#ifdef _FPE_TRAP_
use, intrinsic :: IEEE_Exceptions, only: IEEE_Set_Halting_Mode, IEEE_Usual
use Definitions, only: int32
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: rc

#ifdef _FPE_TRAP_
call IEEE_Set_Halting_Mode(IEEE_Usual,.true._int32)
#endif

call Start('extf')
call extf(rc)
call Finish(rc)

end program Main
