************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2019, Gerardo Raggi                                    *
************************************************************************
      module kriging_mod

      implicit none
!
!     Define and initiate Kriging parameters.
!
      Logical :: Kriging = .False.
      Integer :: nspAI = 1
      Logical :: anMd = .True.
      Real*8  :: pAI = 2
      Real*8  :: lb(3) = [20.0D0, 20.0D0, 1.0D0]
      Integer :: Max_MicroIterations = 50
      Real*8  :: Thr_MicroIterations = 1.0D-8
      Logical :: blAI = .False.
      Logical :: mblAI = .False.
      Logical :: blaAI = .True.
      Real*8  :: blavAI=10.0D0
      Logical :: set_l=.False.
      Logical :: ordinary=.False.
      Real*8  :: blvAI
*
!     Memory for coordinates, value and gradients of the
!     sample points.
!
      real*8, allocatable, protected :: x(:,:), y(:), dy(:)
!
!     Inter_save: the dimension of the coordinate vector
!     nPoints_v: the total number of sample points for which the value is
!                used
!     nPoinst_g: the total number of sample points for which the
!                gradients are used
!
!     We will assume that nPoints_v >= nPoints_g
!
      integer, protected :: nInter_save = 0
      integer, protected :: nPoints_v = 0, nPoints_g = 0

      real*8, allocatable ::
     &        rl(:,:), dl(:), full_Rinv(:,:),
     &        full_R(:,:), x0(:), Kv(:),
     &        cv(:,:,:), cvg(:,:,:),cvh(:,:,:,:),
     &        Rones(:), l(:),
     &        gpred(:), hpred(:,:), ll(:),
     &        cvMatFder(:), cvMatSder(:), cvMatTder(:)
      real*8 :: pred, sigma, var
      real*8 :: sb, variance, detR, lh, sbO, sbmev
      real*8, parameter :: h = 1e-5, eps = 1e-13, eps2 = 1e-10
! eps avoid to become singular in 1st der & eps2 in 2nd der
      integer :: prev_ns, m_t, counttimes

      contains

      Subroutine Setup_Kriging(nPoints_In,nD,nInter,x_,dy_,y_)
#include "stdalloc.fh"
      integer :: nPoints_In, nD, nInter, i, j
      real*8 ::  x_(nInter,nPoints_In)
      real*8 ::         y_(nPoints_In)
      real*8 :: dy_(nInter,nPoints_In-nD)

      nInter_save = nInter

      nPoints_v    = nPoints_In
      nPoints_g    = nPoints_In-nD

      Call mma_Allocate(x,nInter,nPoints_v,Label="x")
      Call mma_Allocate(dy,nInter*nPoints_g,Label="dy")
      Call mma_Allocate(y,nPoints_v,Label="y")

!x is the n-dimensional internal coordinates
      x(:,:) = x_(:,:)
      ! write(6,*) 'x',x
!y is the energy
      y(:) = y_(:)
      ! write(6,*) 'y',y
!dy it's a vector of Grad-y (eq. (5)  ref. gradients of
! the energy with respect to the internal coordinates
!
!     Note the storage as subblocks of the same component of the
!     gradient, each subblock running over all nPoints_g which
!     contributes with gradient values.
!
      do i=1,nInter
        do j=1,nPoints_g
          dy(j + (i-1)*nPoints_g) = dy_(i,j)
        enddo
      enddo

      return
      end subroutine Setup_kriging

      end module kriging_mod
