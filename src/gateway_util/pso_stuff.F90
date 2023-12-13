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

module PSO_Stuff

use Data_structures, only: DSBA_Type
use Definitions, only: wp, iwp

! --- Stuff for Aces 2 read of Gamma file:
! Bin G_Toc Gamma_On lBin LuGamma nGamma nSSDM SO2cI SSDM

implicit none
private

integer(kind=iwp) :: iD0Lbl, iOff_ij2K(8), kCMO, lBin, LuGam, LuGamma, m_Txy, mCMO, mDens, mG1, mG2, n_ij2K, n_Txy, nBasA, &
                     nBasASQ, nBasT, nDens, nG1, nG2, nGamma, nnP(0:7), npos(0:7,3), nSOs1, nSSDM, nV_K, nZ_p_k
logical(kind=iwp) :: lPSO, lsa, Case_3C, Case_2C, Case_mp2, Gamma_mrcisd, Gamma_On, NO_NUC
character(len=7) :: FnGam
integer(kind=iwp), allocatable :: ij2K(:), SO2ci(:,:)
real(kind=wp), allocatable :: A_PT2(:,:), B_PT2(:,:,:), Bin(:,:), CMO(:,:), D0(:,:), DMdiag(:,:), DS(:), DSVar(:), DVar(:,:), &
                              G1(:,:), G2(:,:), G_Toc(:), SSDM(:,:,:), Thpkl(:), Txy(:,:), U_k(:), V_k(:,:), Z_p_k(:,:)
type(DSBA_Type), allocatable, target :: AOrb(:)

public :: A_PT2, AOrb, B_PT2, Bin, Case_2C, Case_3C, Case_mp2, CMO, D0, DMdiag, DS, DSVar, DVar, FnGam, G1, G2, G_Toc, &
          Gamma_mrcisd, Gamma_On, iD0Lbl, ij2K, iOff_ij2K, kCMO, lBin, lPSO, lsa, LuGam, LuGamma, m_Txy, mCMO, mDens, mG1, mG2, &
          n_ij2K, n_Txy, nBasA, nBasASQ, nBasT, nDens, nG1, nG2, nGamma, nnP, NO_NUC, npos, nSOs1, nSSDM, nV_K, nZ_p_k, SO2cI, &
          SSDM, Thpkl, Txy, U_k, V_k, Z_p_k

end module PSO_Stuff
