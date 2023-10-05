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

! Stuff for Cholesky MP2 program (chomp2):
!
! iAOVir, iBatOrb, iDel, iFro, iMatab, iOcc, iT1am, iT1AOT, iVir, lUnit_F, nAOVir, nBatch, nBatOrbT, nDel, nDelT, nFro, nFroT,
! nMatab, nMP2Vec, nOcc, nOccT, nOrb, nPQ_Prod, nT1am, nT1AOT, nTypF, nVir, nVirT, RootNm

! Cholesky MP2 configuration stuff (chomp2_cfg):
!
! all_Vir, C_os, ChkDecoMP2, ChoAlg, Decom_Def, DecoMP2, DoDens, DoFNO, DoGrdt, DoMP2, DoT1amp, EMP2_dens, EOSMP2, FNOMP2,
! ForceBatch, iOffT1, l_Dii, Laplace, Laplace_BlockSize, Laplace_BlockSize_Def, Laplace_nGridPoints, LovMP2, MxQual_Def, MxQualMP2,
! nActa, NoGamma, OED_Thr, set_cd_thr, SOS_mp2, Span_Def, SpanMP2, ThrLov, ThrMP2, Verbose, vkept, Wref, XEMP2
! Laplace_mGridPoints: mx implemented in minimax opt

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

integer(kind=iwp), parameter :: nTypF = 2
integer(kind=iwp) :: ChoAlg, iAdrOff(8,9), iAoMo(8,8,9), iAOVir(8,8), iBatOrb(8), iDel(8), iFro(8), iMatab(8,8), iMoAo(8,8,9), &
                     iMoMo(8,8,9), iOcc(8), iOffT1(8), iOption_MP2CD, iT1am(8,8), iT1AOT(8,8), iVir(8), kFLagr(8), kLagr(8), &
                     kPab(8), kPai(8), kPaK(8), kPij(8), kPiK(8), kWab(8), kWai(8), kWaK(8), kWij(8), kWiK(8), kWJK(8), l_Dii, &
                     Laplace_BlockSize, Laplace_nGridPoints, lFLagr, lLagr, lPab, lPai, lPaK, lPij, lPiK, lUnit_F(8,nTypF), &
                     LuRInv(8), LuUVec, LuVVec, LuWVec, lWab, lWai, lWaK, lWij, lWiK, lWJK, MxQualMP2, nActa, nAdrOff(8), &
                     nAoMo(8,9), nAOVir(8), nBatch, nBatOrbT, nDel(8), nDelT, nFro(8), nFroT, nMatab(8), nMo(8,9), nMoAo(8,9), &
                     nMoMo(8,9), nMoType, nMP2Vec(8), nOcc(8), nOccT, nOccVirT, nOrb(8), NowSym, nPQ_prod(8), nT1am(8), nT1AOT(8), &
                     nVir(8), nVirT
real(kind=wp) :: C_os, DeMP2, EMP2_dens, EOSMP2, OED_Thr, shf, SpanMP2, ThrLov, ThrMP2, vkept, Wref, XEMP2
logical(kind=iwp) :: all_Vir, ChkDecoMP2, ChoMP2_allocated = .false., ChoMP2g_allocated = .false., DecoMP2, DoDens, DoFNO, DoGrdt, &
                     DoMP2, DoT1amp, FNOMP2, ForceBatch, InCore(8), Laplace, LovMP2, MP2_small, NoGamma, set_cd_thr, SOS_mp2, &
                     Verbose
type(V2) :: MP2D(8), MP2D_e(8), MP2W(8), MP2W_e(8)
integer(kind=iwp), allocatable :: AdrR1(:,:,:), AdrR2(:,:,:), iFirst(:), iFirstS(:,:), LiMatij(:,:,:), LiPQprod(:,:,:), &
                                  LiT1am(:,:,:), LnBatOrb(:,:), LnMatij(:,:), LnOcc(:,:), LnPQprod(:,:), LnT1am(:,:), lUnit(:,:), &
                                  NumBatOrb(:), NumOcc(:)
real(kind=wp), allocatable, target :: MP2D_e_full(:), MP2D_full(:), MP2W_e_full(:), MP2W_full(:)
real(kind=wp), allocatable :: EFrozT(:), EOccuT(:), EVirtT(:), OldVec(:), T1amp(:)
real(kind=wp), pointer, contiguous :: EOcc(:) => null(), EVir(:) => null()
integer(kind=iwp), parameter :: Laplace_BlockSize_Def = 500, Laplace_mGridPoints = 20, MxQual_Def = 200
real(kind=wp), parameter :: Span_Def = 1.0e-2_wp
logical(kind=iwp), parameter :: Decom_Def = .false.
character(len=*), parameter :: RootNm = 'ChoMP2_'

public :: AdrR1, AdrR2, all_Vir, C_os, ChkDecoMP2, ChoAlg, ChoMP2_allocated, ChoMP2g_allocated, Decom_Def, DecoMP2, DeMP2, DoDens, &
          DoFNO, DoGrdt, DoMP2, DoT1amp, EFrozT, EMP2_dens, EOcc, EOccuT, EOSMP2, EVir, EVirtT, FNOMP2, ForceBatch, iAdrOff, &
          iAoMo, iAOVir, iBatOrb, iDel, iFirst, iFirstS, iFro, iMatab, iMoAo, iMoMo, InCore, iOcc, iOffT1, iOption_MP2CD, iT1am, &
          iT1AOT, iVir, kFLagr, kLagr, kPab, kPai, kPaK, kPij, kPiK, kWab, kWai, kWaK, kWij, kWiK, kWJK, l_Dii, Laplace, &
          Laplace_BlockSize, Laplace_BlockSize_Def, Laplace_mGridPoints, Laplace_nGridPoints, lFLagr, LiMatij, LiPQprod, LiT1am, &
          lLagr, LnBatOrb, LnMatij, LnOcc, LnPQprod, LnT1am, LovMP2, lPab, lPai, lPaK, lPij, lPiK, lUnit, lUnit_F, LuRInv, LuUVec, &
          LuVVec, LuWVec, lWab, lWai, lWaK, lWij, lWiK, lWJK, MP2_small, MP2D, MP2D_e, MP2D_e_full, MP2D_full, MP2W, MP2W_e, &
          MP2W_e_full, MP2W_full, MxQual_Def, MxQualMP2, nActa, nAdrOff, nAoMo, nAOVir, nBatch, nBatOrbT, nDel, nDelT, nFro, &
          nFroT, nMatab, nMo, nMoAo, nMoMo, nMoType, NMP2Vec, nOcc, nOccT, nOccVirT, NoGamma, nOrb, NowSym, nPQ_Prod, nT1am, &
          nT1AOT, nTypF, NumBatOrb, NumOcc, nVir, nVirT, OED_Thr, OldVec, RootNm, Set_cd_thr, shf, SOS_mp2, Span_Def, SpanMP2, &
          T1amp, ThrLov, ThrMP2, Verbose, vkept, Wref, XEMP2

end module ChoMP2
