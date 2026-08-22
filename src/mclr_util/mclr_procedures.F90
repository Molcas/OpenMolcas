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
module MCLR_procedures

implicit none
private

public :: CISigma, CISigma_sa, CISigma_td, RdJobIph, RdJobIph_td

contains

#define _IN_MODULE_
#include "../mclr/cisigma.F90"
#include "../mclr/cisigma_sa.F90"
#include "../mclr/cisigma_td.F90"
#include "../mclr/rdjobiph.F90"
#include "../mclr/rdjobiph_td.F90"

end module MCLR_procedures
