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
Module nq_Grid
Real*8, Allocatable:: Weights(:)
Real*8, Allocatable:: Grid(:,:)
!     nGridMax: size of the array Grid
Integer :: nGridMax=128
Real*8, Allocatable:: Rho(:,:)
Real*8, Allocatable:: vRho(:,:)
Integer :: nRho=0
Real*8, Allocatable:: GradRho(:,:)
Integer :: nGradRho=0
Real*8, Allocatable:: Sigma(:,:)
Real*8, Allocatable:: vSigma(:,:)
Integer :: nSigma=0
Real*8, Allocatable:: Lapl(:,:)
Real*8, Allocatable:: vLapl(:,:)
Integer :: nLapl=0
Real*8, Allocatable:: Tau(:,:)
Real*8, Allocatable:: vTau(:,:)
Integer :: nTau=0
Logical :: l_CASDFT=.FALSE.
Real*8, Allocatable:: Exc(:)
Real*8, Allocatable:: F_xca(:), F_xcb(:)
Real*8, Allocatable, Target:: TabAO(:,:,:)
Real*8, Allocatable, Target:: TabAO_Short(:,:,:)
Real*8, Pointer:: TabAO_pack(:) => Null()
Real*8, Allocatable:: Grid_AO(:,:,:,:)
Real*8, Allocatable:: Dens_AO(:,:,:)
Real*8, Allocatable:: dRho_dR(:,:,:)
Integer, Allocatable:: Ind_Grd(:,:)
Integer, Allocatable:: iBfn_Index(:,:)
Integer :: kAO=0
End Module nq_Grid
