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

module DMRGInfo
! These are used for the MCLR part
! ndets_RGLR   : number of Slater determinants
! LRras2       : CI space when solving LR equation
! RGras2       : CI space when solving LR equation
! nstates_RGLR : number of states

use Definitions, only: iwp

implicit none
private

logical(kind=iwp) :: doDMRG, doMCLR

integer(kind=iwp) :: LRras2(8), MS2_RGLR, ndets_RGLR, nele_RGLR, nstates_RGLR, RGras2(8)

public :: DoDMRG, DoMCLR, LRras2, MS2_RGLR, ndets_RGLR, nele_RGLR, nstates_RGLR, RGras2

end module DMRGInfo
