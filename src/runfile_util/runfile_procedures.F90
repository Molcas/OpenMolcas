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
module RunFile_procedures

implicit none
private

public :: Get_Coord_New, Get_dExcdRa, Get_PC_Coord_New

contains

#define _IN_MODULE_
#include "get_coord_new.F90"
#include "get_pc_coord_new.F90"
#include "get_dexcdra.F90"

end module RunFile_procedures
