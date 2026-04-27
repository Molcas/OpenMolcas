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

subroutine Initialize_BLAS(lib,prlev)

use link_blas, only: lb_initialize
use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: lib
integer(kind=iwp), intent(in) :: prlev

call lb_initialize(lib,prlev)

end subroutine Initialize_BLAS
