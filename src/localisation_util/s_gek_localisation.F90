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
! Copyright (C) 2026, Lila Zapp                                        *
!***********************************************************************

subroutine S_GEK_localisation()

use Definitions, only: wp, iwp, u6

implicit none

integer(kind=iwp) :: mDiis, nDiis, Max_Iter
real(kind=wp), intent(inout) :: q_diis(mDiis,nDiis+Max_Iter), g_diis(mDiis,nDiis+Max_Iter), Energy(nDiis+Max_Iter), &
                                H_diis(mDiis,mDiis)
real(kind=wp) :: dq_diis(mDiis), dqdq
character, intent(out) :: Step_Trunc
character(len=6), intent(out) :: UpMeth
logical(kind=iwp), intent(in) :: SORange

integer(kind=iwp), intent(in) :: mOV
real(kind=wp), intent(inout) :: dq(mOV)
real(kind=wp), intent(out) :: dqdq
character(len=6), intent(inout) :: UpMeth
character, intent(inout) :: Step_Trunc
logical(kind=iwp), intent(in) :: SOrange
integer(kind=iwp) :: i, iFirst, ipg, ipq, j, k, l, mDIIS, nDIIS, nExplicit
real(kind=wp) :: Cpu1, Cpu2, gg, Tim1, Tim2, Tim3
real(kind=wp), allocatable :: D(:,:), dq_diis(:), e_diis(:,:), g(:,:), g_diis(:,:), H_Diis(:,:), q(:,:), q_diis(:,:), w(:,:)
integer(kind=iwp), parameter :: Max_Iter = 50, nWindow = 20
real(kind=wp), parameter :: Beta_Disp_Min = 5.0e-3_wp, Beta_Disp_Seed = 0.05_wp, StepMax_Seed = 0.1_wp, Thr_RS = 1.0e-7_wp, &
                            ThrGrd = 1.0e-7_wp
#ifndef _FULL_SPACE_
real(kind=wp), allocatable :: aux_a(:), aux_b(:)
#endif
real(kind=wp), external :: DDot_


! new stuff added for the GEK
integer(kind=iwp) ::fullspace_dim,&
                    n_SGEK, & !n^diis ()
                    m_SGEK, &
                    nExplicit ! number of vectors used to get a basis (some are likely linearly dependent)
integer(kind=iwp), parameter :: Max_Iter_SGEK = 50, nWindow_SGEK = 20



!subroutine GEK_Optimizer(mDiis,nDiis,Max_Iter,q_diis,g_diis,dq_diis,Energy,H_diis,dqdq,Step_Trunc,UpMeth,SORange)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! mDiis             subspace dimensionality (<=2*ndiis); number of linear independent e_diis column vectors
! nDiis             number of iterations used to span the subspace; nDIIS = min(IterGEK,nWindow)
! Max_Iter          maximal number of iterations to consider for the GEK surrogate model; usually Max_Iter = 50
! q_diis            projected coordinate vectors
! g_diis            projected gradient vectors
! dq_diis           output displacement, suggested by the optimization in the subspace; still in subspace representation
!                   get fullspace representation by doing: dq(:) = dq(:)+dq_diis(i)*e_diis(:,i) for i=1,mdiis
! Energy            y vector
! H_diis            projected Hessian diagonal
! dqdq              some output, (real) that has to do with the full space displacement
! Step_Trunc        some output (character)
! UpMeth            some output (string), e.g. "RVO"
! SORange           some input (logical)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the last n_SGEK iterations
!n_SGEK = min(nWindow_SGEK,Iter) ! e.g. first iteration Iter =1, m_SGEK = 1

!call GEK_Optimizer(mDiis,nDiis,Max_Iter,q_diis,g_diis,dq_diis,Energy,H_diis,dqdq,Step_Trunc,UpMeth,SORange)


end subroutine S_GEK_localisation
