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

! This module contains procedures that need an interface
module Do_Grid

implicit none
private

public :: Do_GGL, Do_Lebedev, Do_Lebedev_Sym, Do_Lobatto

contains

#define _IN_MODULE_
#include "do_ggl.F90"
#include "do_lebedev.F90"
#include "do_lebedev_sym.F90"
#include "do_lobatto.F90"

end module Do_Grid
