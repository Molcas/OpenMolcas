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

! Wrapper module to allow using the subroutines without explicit interface

module Fock_util_interface

implicit none
private

public :: cho_lr_mos, choras_drv

contains

#include "cho_lr_mos.F90"
#include "choras_drv.F90"

end module Fock_util_interface
