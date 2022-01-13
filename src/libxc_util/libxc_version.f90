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
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************
Subroutine libxc_version()
  use xc_f03_lib_m
  implicit none
  integer*4 :: vmajor, vminor, vmicro
  ! Print out the version
  call xc_f03_version(vmajor, vminor, vmicro)
  write(6,'(6X,"Libxc version: ",I1,".",I1,".",I1)') vmajor, vminor, vmicro
End Subroutine libxc_version
