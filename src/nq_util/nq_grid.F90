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

module nq_Grid

use Definitions, only: wp, iwp

implicit none
private

! nGridMax: size of the array Grid
integer(kind=iwp) :: kAO = 0, nGridMax = 128, nRho = 0
logical(kind=iwp) :: l_CASDFT = .false.
integer(kind=iwp), allocatable :: Angular(:), iBfn_Index(:,:), IndGrd(:), iTab(:,:), List_G(:,:), nR_Eff(:)
real(kind=wp), allocatable :: Coor(:,:), D1Unzip(:,:), Dens_AO(:,:,:), dRho_dR(:,:,:), dW_dR(:,:), F_xc(:), F_xca(:), F_xcb(:), &
                              Fact(:,:), GradRho(:,:), Grid(:,:), Grid_AO(:,:,:,:), Lapl(:,:), Mem(:), P2Unzip(:,:,:,:), Pax(:,:), &
                              R2_trial(:), Rho(:,:), Sigma(:,:), Tau(:,:), Temp(:), vLapl(:,:), vRho(:,:), vSigma(:,:), vTau(:,:), &
                              Weights(:)
real(kind=wp), allocatable, target :: TabAO(:,:,:), TabAO_Short(:,:,:)
real(kind=wp), pointer :: TabAO_pack(:) => null()

public :: Angular, Coor, D1Unzip, Dens_AO, dRho_dR, dW_dR, F_xc, F_xca, F_xcb, Fact, GradRho, Grid, Grid_AO, iBfn_Index, IndGrd, &
          iTab, kAO, l_CASDFT, Lapl, List_G, Mem, nGridMax, nR_Eff, nRho, P2Unzip, Pax, R2_trial, Rho, Sigma, TabAO, TabAO_pack, &
          TabAO_Short, Tau, Temp, vLapl, vRho, vSigma, vTau, Weights

end module nq_Grid
