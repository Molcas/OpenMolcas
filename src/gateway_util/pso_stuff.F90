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

logical lPSO, lsa, Case_3C, Case_2C, Case_mp2

integer nnP(0:7), iOff_ij2K(8), npos(0:7,3)

real*8, allocatable :: DMdiag(:,:), Thpkl(:), G2(:,:), CMO(:,:)
real*8, allocatable :: Txy(:,:), V_k(:,:), U_k(:), Z_p_k(:,:)
real*8, allocatable :: G1(:,:), D0(:,:), DVar(:,:), DS(:), DSVar(:)
integer, allocatable :: ij2K(:)

integer nG2, mG2
integer nG1, mG1
integer mCMO, kCMO
integer nDens, mDens
integer n_Txy, m_Txy
integer n_ij2K
integer nZ_p_k
integer nV_K, nSOs1
integer iD0Lbl

type(DSBA_Type), allocatable, target :: AOrb(:)

end module PSO_Stuff
