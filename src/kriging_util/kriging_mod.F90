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
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************

module kriging_mod

use Definitions, only: wp, iwp

implicit none
private

! Define and initiate Kriging parameters.

logical(kind=iwp) :: Kriging = .false., &
                     anMd = .true., &
                     blAI = .false., &
                     mblAI = .false., &
                     blaAI = .true., &
                     set_l = .false., &
                     ordinary = .false., &
                     PGEK_On = .false.
integer(kind=iwp) :: nspAI = 1, pAI = 2, Max_MicroIterations = 50, nD_In = 0
real(kind=wp) :: lb(3) = [20.0_wp,20.0_wp,1.0_wp], &
                 Thr_MicroIterations = 1.0e-8_wp, &
                 blavAI = 10.0_wp, &
                 blvAI = 0.0_wp

! Memory for coordinates, value and gradients of the
! sample points.

real(kind=wp), allocatable, protected :: x(:,:), y(:), dy(:)

! Inter  : the dimension of the coordinate vector
! nPoints: the total number of sample points for which the value is
!          used or the gradient is used
! nD     : the total number of sample points less for which the
!          gradients are used
!
! We will assume that nD >= 0

integer(kind=iwp), protected :: nInter = 0, nPoints = 0, nD = 0
integer(kind=iwp) :: nInter_Eff = 0

real(kind=wp), allocatable :: rl(:,:), dl(:), full_Rinv(:,:), full_R(:,:), x0(:), Kv(:), cv(:,:,:), Rones(:), l(:), gpred(:), &
                              hpred(:,:), ll(:), cvMatFder(:), cvMatSder(:), cvMatTder(:)
integer(kind=iwp), allocatable :: Index_PGEK(:)
real(kind=wp) :: pred, sigma, var, sb, variance, detR, lh, sbO, sbmev
integer(kind=iwp) :: m_t
real(kind=wp), parameter :: h = 1e-5, eps = 1e-13, eps2 = 1e-10
! eps avoid to become singular in 1st der & eps2 in 2nd der

! Transformation matrix for kriging_layer routines

real(kind=wp), allocatable :: layer_U(:,:)

public :: anMd, blaAI, blAI, blavAI, blvAI, cv, cvMatFder, cvMatSder, cvMatTder, detR, dl, dy, eps, eps2, full_R, full_RInv, &
          gpred, h, hpred, Index_PGEK, Kriging, kv, l, layer_U, lb, lh, ll, Max_MicroIterations, m_t, mblAI, nD, nD_In,  nInter, &
          nInter_Eff, nPoints, nspAI, ordinary, pAI, PGEK_On, pred, rl, Rones, sb, sbmev, sbO, set_l, sigma, Thr_MicroIterations, &
          var, variance, x, x0, y
public :: Deallocate_Protected, Prep_Kriging

contains

subroutine Prep_Kriging(nPoints_In,nInter_In,x_,dy_,y_)

  use stdalloc, only: mma_allocate

  integer(kind=iwp), intent(in) :: nPoints_In, nInter_In
  real(kind=wp), intent(in) :: x_(nInter_In,nPoints_In), y_(nPoints_In), dy_(nInter_In,nPoints_In)
  integer(kind=iwp) :: i, j

  nInter = nInter_In
  nInter_Eff = nInter
  nPoints = nPoints_In
  nD = max(0,min(nD_In,nPoints-nD_In))

  ! Allocate arrays for data or energies, coordinates, and gradients

  call mma_Allocate(x,nInter,nPoints,label='x')
  call mma_Allocate(y,nPoints,label='y')
  call mma_Allocate(dy,nInter*(nPoints-nD),label='dy')

  ! The code will use partial GEK with indirect addressing. However,
  ! here we defaults the index array so that it behaves as conventional GEK.

  call mma_Allocate(Index_PGEK,nInter,label='Index_PGEK')
  do i=1,nInter
    Index_PGEK(i) = i
  end do

  ! x is the n-dimensional internal coordinates
  x(:,:) = x_(:,:)
  !write(u6,*) 'x',x
  ! y is the energy
  y(:) = y_(:)
  !write(u6,*) 'y',y
  ! dy is a vector of Grad-y (eq. (5) ref. gradients of
  ! the energy with respect to the internal coordinates
  !
  ! Note the storage as subblocks of the same component of the
  ! gradient, each subblock running over all nPoints_g which
  ! contributes with gradient values.
  ! This will enable a somewhat simpler code later on in the PGEK case.
  !
  ! At this point we also skip those gradients which we will not use in
  ! the kriging.

  do i=1,nInter
    do j=1,nPoints-nD
      dy(j+(i-1)*(nPoints-nD)) = dy_(i,j+nD)
    end do
  end do

  return

end subroutine Prep_kriging

subroutine Deallocate_protected()

  use stdalloc, only: mma_deallocate

  call mma_Deallocate(x)
  call mma_Deallocate(y)
  call mma_Deallocate(dy)

end subroutine Deallocate_protected

end module kriging_mod
