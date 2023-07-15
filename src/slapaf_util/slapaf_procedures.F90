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
module Slapaf_procedures

private
public :: Box, Hidden, OldFCM, SphInt

contains

#define _IN_MODULE_
#include "box.F90"
#include "hidden.F90"
#include "oldfcm.F90"

end module Slapaf_procedures
