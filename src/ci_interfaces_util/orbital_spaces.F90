!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

!> The orbital spaces per symmetry (frozen, inactive, active, secondary,
!> deleted; basis functions and orbitals) and their totals.
module Orbital_Spaces

use Molcas, only: MxSym
use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: NASH(mxSym), NBAS(mxSym), NDEL(mxSym), NDELT, NFRO(mxSym), NFROT, NISH(mxSym), NORB(mxSym), NSSH(mxSym), &
                     NSYM, NTOT, NTOT1, NTOT2

public :: NASH, NBAS, NDEL, NDELT, NFRO, NFROT, NISH, NORB, NSSH, NSYM, NTOT, NTOT1, NTOT2

end module Orbital_Spaces
