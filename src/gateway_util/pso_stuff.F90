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

implicit none
private

integer(kind=iwp) :: iD0Lbl, iOff_ij2K(8), kCMO, m_Txy, mCMO, mDens, mG1, mG2, n_ij2K, n_Txy, nBasA, nBasASQ, nBasT, nDens, nG1, &
                     nG2, nnP(0:7), npos(0:7,3), nSOs1, nV_K, nZ_p_k
logical(kind=iwp) :: lPSO, lsa, Case_3C, Case_2C, Case_mp2
integer(kind=iwp), allocatable :: ij2K(:)
real(kind=wp), allocatable :: CMO(:,:), D0(:,:), DMdiag(:,:), DS(:), DSVar(:), DVar(:,:), G1(:,:), G2(:,:), Thpkl(:), Txy(:,:), &
                              U_k(:), V_k(:,:), Z_p_k(:,:), A_PT2(:,:), B_PT2(:,:,:)
type(DSBA_Type), allocatable, target :: AOrb(:)

public :: A_PT2, AOrb, B_PT2, Case_2C, Case_3C, Case_mp2, CMO, D0, DMdiag, DS, DSVar, DVar, G1, G2, iD0Lbl, ij2K, iOff_ij2K, kCMO, &
          lPSO, lsa, m_Txy, mCMO, mDens, mG1, mG2, n_ij2K, n_Txy, nBasA, nBasASQ, nBasT, nDens, nG1, nG2, nnP, npos, nSOs1, nV_K, &
          nZ_p_k, Thpkl, Txy, U_k, V_k, Z_p_k

end module PSO_Stuff
