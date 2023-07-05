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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************

module ThermoChem

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: nsRot = 0, nUserPT = 0
real(kind=wp) :: UserP = Zero, UserT(64) = Zero
logical(kind=iwp) :: lDoubleIso = .false., lTherm = .false.

public :: lDoubleIso, lTherm, nsRot, nUserPT, UserP, UserT

end module ThermoChem
