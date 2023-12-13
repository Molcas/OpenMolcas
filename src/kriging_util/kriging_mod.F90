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

use Constants, only: Zero
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
real(kind=wp) :: Thr_MicroIterations = 1.0e-8_wp, &
                 blavAI = 10.0_wp, &
                 blvAI = Zero

! Memory for coordinates, value and gradients of the
! sample points.

real(kind=wp), allocatable, protected :: x(:,:), y(:,:), dy(:,:)

! Inter  : the dimension of the coordinate vector
! nPoints: the total number of sample points for which the value is
!          used or the gradient is used
! nD     : the total number of sample points less for which the
!          gradients are used
!
! We will assume that nD >= 0

integer(kind=iwp), protected :: nInter = 0, nPoints = 0, nD = 0
integer(kind=iwp) :: nInter_Eff = 0, nSet = 1

integer(kind=iwp) :: m_t
real(kind=wp) :: detR, sbmev, sbO, var
integer(kind=iwp), allocatable :: Index_PGEK(:), Model_Type(:)
real(kind=wp), allocatable :: cv(:,:,:), cvMatFder(:), cvMatSder(:), cvMatTder(:), dl(:), full_R(:,:), full_Rinv(:,:), gpred(:,:), &
                              hpred(:,:,:), Kv(:,:), l(:), lh(:), pred(:), rl(:,:), Rones(:), sb(:), sigma(:), variance(:), x0(:)
! eps avoid to become singular in 1st der & eps2 in 2nd der
real(kind=wp), parameter :: h = 1.0e-5_wp, eps = 1.0e-13_wp, eps2 = 1.0e-10_wp

! Transformation matrix for kriging_layer routines
real(kind=wp), allocatable :: layer_U(:,:)

public :: anMd, blaAI, blAI, blavAI, blvAI, cv, cvMatFder, cvMatSder, cvMatTder, detR, dl, dy, eps, eps2, full_R, full_RInv, &
          gpred, h, hpred, Index_PGEK, Kriging, kv, l, layer_U, lh, Max_MicroIterations, m_t, mblAI, Model_Type, nD, nD_In, &
          nInter, nInter_Eff, nPoints, nSet, nspAI, ordinary, pAI, PGEK_On, pred, rl, Rones, sb, sbmev, sbO, set_l, sigma, &
          Thr_MicroIterations, var, variance, x, x0, y
public :: Deallocate_Protected, Prep_Kriging

contains

subroutine Prep_Kriging(nPoints_In,nInter_In,x_,dy_,y_)

  use stdalloc, only: mma_allocate

  integer(kind=iwp), intent(in) :: nPoints_In, nInter_In
  real(kind=wp), intent(in) :: x_(nInter_In,nPoints_In), y_(nPoints_In,nSet), dy_(nInter_In,nPoints_In,nSet)
  integer(kind=iwp) :: i, j

  nInter = nInter_In
  nInter_Eff = nInter
  nPoints = nPoints_In
  nD = max(0,min(nD_In,nPoints-nD_In))

  ! Allocate arrays for data or energies, coordinates, and gradients

  call mma_Allocate(x,nInter,nPoints,label='x')
  call mma_Allocate(y,nPoints,nSet,label='y')
  call mma_Allocate(dy,nInter*(nPoints-nD),nSet,label='dy')

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
  y(:,:) = y_(:,:)
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
      dy(j+(i-1)*(nPoints-nD),:) = dy_(i,j+nD,:)
    end do
  end do

  return

end subroutine Prep_kriging

subroutine Deallocate_protected()

  use stdalloc, only: mma_deallocate

  call mma_Deallocate(x)
  call mma_Deallocate(y)
  call mma_Deallocate(dy)
  if (allocated(Model_Type)) call mma_Deallocate(Model_Type)

end subroutine Deallocate_protected

end module kriging_mod
