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


module Surfacehop_globals

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: InitSeed, NSUBSTEPS
real(kind=wp) :: DECO, Ethreshold, FixedRand, RandThreshold
logical(kind=iwp) :: decoherence, firststep, fixedrandL, iseedL, lH5Restart, rassi_ovlp, Run_rassi, tullyL, tullySubVerb
character(len=180) :: File_H5Res

public :: DECO, decoherence, Ethreshold, File_H5Res, firststep, FixedRand, fixedrandL, InitSeed, iseedL, lH5Restart, NSUBSTEPS, &
          RandThreshold, rassi_ovlp, Run_rassi, tullyL, tullySubVerb

end module Surfacehop_globals
