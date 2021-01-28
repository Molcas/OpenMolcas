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

module Tully_variables

use Definitions, only: wp, iwp

implicit none
private

logical(kind=iwp) :: tullyL, decoherence, tullySubVerb, fixedrandL, iseedL
real(kind=wp) :: DECO, Ethreshold, RandThreshold, FixedRand
integer(kind=iwp) :: NSUBSTEPS, InitSeed

public :: tullyL, decoherence, tullySubVerb, fixedrandL, iseedL, DECO, Ethreshold, RandThreshold, FixedRand, NSUBSTEPS, InitSeed

end module
