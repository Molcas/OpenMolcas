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
module RI_procedures

implicit none
private

public :: Drv2El_2Center_RI, Drv2el_Atomic_Nosym, Effective_CD_Pairs, Fix_Exponents

contains

#define _IN_MODULE_
#include "drv2el_2center_ri.F90"
#include "drv2el_atomic_nosym.F90"
#include "effective_cd_pairs.F90"
#include "fix_exponents.F90"

end module RI_procedures
