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

program main

#ifdef _FPE_TRAP_
use, intrinsic :: ieee_exceptions
#endif
implicit real*8(a-h,o-z)
character*20 Module_Name
parameter(Module_Name='localisation')

#ifdef _FPE_TRAP_
call IEEE_Set_Halting_Mode(IEEE_Usual,.true._4)
#endif

call Start(Module_Name)
call Localisation(ireturn)
call Finish(ireturn)

end program main
