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
  use, intrinsic :: IEEE_Exceptions
#endif
  implicit none
  integer :: ireturn
  character(len=*), parameter :: module_name = 'rhodyn'
#ifdef _FPE_TRAP_
  call IEEE_Set_Halting_Mode(IEEE_Usual,.true._4)
#endif

  call start(module_name)
  call rhodyn(ireturn)
  call finish(ireturn)
end program main
