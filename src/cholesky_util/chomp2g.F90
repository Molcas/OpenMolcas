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

public :: nMoMo, iMoMo, nMoType
public :: iMoAo, nMoAo, iAoMo, nAoMo
public :: iAdrOff, nAdrOff
public :: nOccVirT
public :: kPab, kPij, kPiK, kLagr, kFLagr, kWab, kWij, kWai, kWiK, kWaK, kPaK, kPai, kWJK, lPab, lPij, lPiK, lLagr, lFLagr, lWab, &
          lWij, lWai, lWiK, lWaK, lPaK, lPai, lWJK
public :: LuUVec, LuVVec, LuWVec, LuRInv

integer nMoMo(8,9), nMo(8,9), iMoMo(8,8,9), nMoType, nMoAo(8,9), iMoAo(8,8,9), nAoMo(8,9), iAoMo(8,8,9), iAdrOff(8,9), nAdrOff(8), &
        nOccVirT, kPab(8), kPij(8), kPiK(8), kLagr(8), kFLagr(8), kWab(8), kWij(8), kWai(8), kWiK(8), kWaK(8), kPaK(8), kPai(8), &
        kWJK(8), lPab, lPij, lPiK, lLagr, lWab, lWij, lWiK, lWaK, lWai, lPaK, lPai, lWJK, lFLagr, LuUVec, LuVVec, LuWVec, LuRInv(8)

end module ChoMP2g
