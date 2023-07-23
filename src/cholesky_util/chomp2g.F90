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

! Common space for Cholesky Mp2-gradients
module ChoMP2g

use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: iAdrOff(8,9), iAoMo(8,8,9), iMoAo(8,8,9), iMoMo(8,8,9), kFLagr(8), kLagr(8), kPab(8), kPai(8), kPaK(8), &
                     kPij(8), kPiK(8), kWab(8), kWai(8), kWaK(8), kWij(8), kWiK(8), kWJK(8), lFLagr, lLagr, lPab, lPai, lPaK, &
                     lPij, lPiK, LuRInv(8), LuUVec, LuVVec, LuWVec, lWab, lWai, lWaK, lWij, lWiK, lWJK, nAdrOff(8), nAoMo(8,9), &
                     nMo(8,9), nMoAo(8,9), nMoMo(8,9), nMoType, nOccVirT

public :: iAdrOff, iAoMo, iMoAo, iMoMo, kFLagr, kLagr, kPab, kPai, kPaK, kPij, kPiK, kWab, kWai, kWaK, kWij, kWiK, kWJK, lFLagr, &
          lLagr, lPab, lPai, lPaK, lPij, lPiK, LuRInv, LuUVec, LuVVec, LuWVec, lWab, lWai, lWaK, lWij, lWiK, lWJK, nAdrOff, nAoMo, &
          nMo, nMoAo, nMoMo, nMoType, nOccVirT

end module ChoMP2g
