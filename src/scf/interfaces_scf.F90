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

module Interfaces_SCF

implicit none
private

public :: dOne_SCF, MinDns, OccDef, PMat_SCF, TraClc_i, vOO2OV

contains

! Subroutines that need an explicit interface (target or allocatable arguments)
#define _IN_MODULE_
#include "done_scf.F90"
#include "mindns.F90"
#include "occdef.F90"
#include "pmat_scf.F90"
#include "traclc_i.F90"
#include "voo2ov.F90"

end module Interfaces_SCF
