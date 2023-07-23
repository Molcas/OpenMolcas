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

! Stuff for decomposing (ai|bj) integrals or amplitudes in MP2:
module chomp2_dec

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: iOption_MP2CD, NowSym
real(kind=wp), pointer, contiguous :: EOcc(:) => null(), EVir(:) => null()
logical(kind=iwp) :: InCore(8)

public :: EOcc, EVir, InCore, iOption_MP2CD, NowSym

end module chomp2_dec
