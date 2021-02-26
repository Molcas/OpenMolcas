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

module info_expbas_mod

use Definitions, only: iwp

implicit none
private

!> Different kinds of orbitals
!> f, i, 1, 2, 3, s, d
integer(kind=iwp), parameter :: n_orb_kinds = 7

logical(kind=iwp) :: DoExpbas, DoDesy
character(len=512) :: EB_FileOrb

public :: DoExpbas, DoDesy, EB_FileOrb, n_orb_kinds

! Exporting some useful parameters
#include "Molcas.fh"
public :: LenIn, MxAtom, mxsym

end module info_expbas_mod
