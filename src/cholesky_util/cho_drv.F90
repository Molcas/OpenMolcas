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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_Drv(iReturn)
!
! Thomas Bondo Pedersen, April 2010.
!
! Purpose: driver for the Cholesky decomposition of two-electron
!          integrals. On entry, the integral program must have been
!          initialized and the relevant index info (#irreps, basis
!          functions, shells, etc.) must have been set up.

use Cholesky, only: Cho_DecAlg
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: iReturn

iReturn = 0
if (Cho_DecAlg == 5) then ! parallel two-step algorithm
  call Cho_Drv_ParTwoStep(iReturn)
else
  call Cho_Drv_Inner(iReturn)
end if

end subroutine Cho_Drv
