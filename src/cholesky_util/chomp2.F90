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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************

module ChoMP2

use Data_Structures, only: V2
use Definitions, only: wp, iwp

! for Cholesky Mp2-gradients (chomp2g):
!
! iAdrOff, iAoMo, iMoAo, iMoMo, kFLagr, kLagr, kPab, kPai, kPaK, kPij, kPiK, kWab, kWai, kWaK, kWij, kWiK, kWJK, lFLagr, lLagr,
! lPab, lPai, lPaK, lPij, lPiK, LuRInv, LuUVec, LuVVec, LuWVec, lWab, lWai, lWaK, lWij, lWiK, lWJK, nAdrOff, nAoMo, nMo, nMoAo,
! nMoMo, nMoType, nOccVirT

! Stuff for decomposing (ai|bj) integrals or amplitudes in MP2 (chomp2_dec):
!
! EOcc, EVir, InCore, iOption_MP2CD, NowSym

implicit none
private

integer(kind=iwp) :: iAdrOff(8,9), iAoMo(8,8,9), iMoAo(8,8,9), iMoMo(8,8,9), iOption_MP2CD, kFLagr(8), kLagr(8), kPab(8), kPai(8), &
                     kPaK(8), kPij(8), kPiK(8), kWab(8), kWai(8), kWaK(8), kWij(8), kWiK(8), kWJK(8), lFLagr, lLagr, lPab, lPai, &
                     lPaK, lPij, lPiK, LuRInv(8), LuUVec, LuVVec, LuWVec, lWab, lWai, lWaK, lWij, lWiK, lWJK, nAdrOff(8), &
                     nAoMo(8,9), nMo(8,9), nMoAo(8,9), nMoMo(8,9), nMoType, nOccVirT, NowSym
logical(kind=iwp) :: ChoMP2_allocated = .false., ChoMP2g_allocated = .false., InCore(8)
type(V2) :: MP2D(8), MP2D_e(8), MP2W(8), MP2W_e(8)
integer(kind=iwp), allocatable :: AdrR1(:,:,:), AdrR2(:,:,:), iFirst(:), iFirstS(:,:), LiMatij(:,:,:), LiPQprod(:,:,:), &
                                  LiT1am(:,:,:), LnBatOrb(:,:), LnMatij(:,:), LnOcc(:,:), LnPQprod(:,:), LnT1am(:,:), lUnit(:,:), &
                                  NumBatOrb(:), NumOcc(:)
real(kind=wp), allocatable, target :: MP2D_e_full(:), MP2D_full(:), MP2W_e_full(:), MP2W_full(:)
real(kind=wp), allocatable :: EFrozT(:), EOccuT(:), EVirtT(:), OldVec(:)
real(kind=wp), pointer, contiguous :: EOcc(:) => null(), EVir(:) => null()

public :: AdrR1, AdrR2, ChoMP2_allocated, ChoMP2g_allocated, EFrozT, EOcc, EOccuT, EVir, EVirtT, iAdrOff, iAoMo, iFirst, iFirstS, &
          iMoAo, iMoMo, InCore, iOption_MP2CD, kFLagr, kLagr, kPab, kPai, kPaK, kPij, kPiK, kWab, kWai, kWaK, kWij, kWiK, kWJK, &
          lFLagr, LiMatij, LiPQprod, LiT1am, lLagr, LnBatOrb, LnMatij, LnOcc, LnPQprod, LnT1am, lPab, lPai, lPaK, lPij, lPiK, &
          lUnit, LuRInv, LuUVec, LuVVec, LuWVec, lWab, lWai, lWaK, lWij, lWiK, lWJK, MP2D, MP2D_e, MP2D_e_full, MP2D_full, MP2W, &
          MP2W_e, MP2W_e_full, MP2W_full, nAdrOff, nAoMo, nMo, nMoAo, nMoMo, nMoType, nOccVirT, NowSym, NumBatOrb, NumOcc, OldVec

end module ChoMP2
