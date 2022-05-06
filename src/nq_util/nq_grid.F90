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

real*8, allocatable :: Pax(:,:)
real*8, allocatable :: Coor(:,:)
real*8, allocatable :: R2_trial(:)
real*8, allocatable :: Fact(:,:)
real*8, allocatable :: Mem(:)
integer, allocatable :: Angular(:)
integer, allocatable :: nR_Eff(:)

integer, allocatable :: List_G(:,:)
integer, allocatable :: iTab(:,:)
integer, allocatable :: IndGrd(:)
real*8, allocatable :: Temp(:)
real*8, allocatable :: P2Unzip(:), D1Unzip(:)
real*8, allocatable :: dW_dR(:,:)

real*8, allocatable :: Weights(:)
real*8, allocatable :: Grid(:,:)
! nGridMax: size of the array Grid
integer :: nGridMax = 128
real*8, allocatable :: Rho(:,:)
real*8, allocatable :: vRho(:,:)
integer :: nRho = 0
real*8, allocatable :: GradRho(:,:)
integer :: nGradRho = 0
real*8, allocatable :: Sigma(:,:)
real*8, allocatable :: vSigma(:,:)
integer :: nSigma = 0
real*8, allocatable :: Lapl(:,:)
real*8, allocatable :: vLapl(:,:)
integer :: nLapl = 0
real*8, allocatable :: Tau(:,:)
real*8, allocatable :: vTau(:,:)
integer :: nTau = 0
logical :: l_CASDFT = .false.
real*8, allocatable :: F_xc(:), F_xca(:), F_xcb(:)
real*8, allocatable, target :: TabAO(:,:,:)
real*8, allocatable, target :: TabAO_Short(:,:,:)
real*8, pointer :: TabAO_pack(:) => null()
real*8, allocatable :: Grid_AO(:,:,:,:)
real*8, allocatable :: Dens_AO(:,:,:)
real*8, allocatable :: dRho_dR(:,:,:)
integer, allocatable :: iBfn_Index(:,:)
integer :: kAO = 0

end module nq_Grid
