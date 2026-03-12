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

subroutine initialize_blas(lib,prlev)
  use link_blas
  implicit none
  character(len=*), intent(in) :: lib
  integer, intent(in) :: prlev
  call lb_initialize(lib,prlev)
end subroutine initialize_blas

subroutine close_blas()
  use link_blas
  implicit none
  call lb_close()
end subroutine close_blas
