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
      Integer :: npxAI = 1
      Real*8  :: lb(3) = [20.0D0, 20.0D0, 1.0D0]
      Integer :: miAI = 50
      Real*8  :: meAI = 1.0D-8
      Logical :: blAI = .False.
      Logical :: mblAI = .False.
      Logical :: blaAI = .True.
      Real*8  :: blavAI=10.0D0
      Logical :: set_l=.False.
      Logical :: ordinary=.False.
      Real*8  :: blvAI
*
      real*8, allocatable, protected :: x(:,:), y(:), dy(:)
      integer, protected :: nInter_save = 0, nPoints_save = 0

      real*8, allocatable ::
     &        rl(:,:,:), dl(:,:), full_Rinv(:,:),
     &        full_R(:,:), nx(:,:), Kv(:),
     &        cv(:,:,:,:), cvg(:,:,:),cvh(:,:,:,:),
     &        var(:), Rones(:), sigma(:), l(:),
     &        pred(:), gpred(:,:), hpred(:,:,:), ll(:),
     &        cvMatFder(:,:), cvMatSder(:,:), cvMatTder(:,:)
      real*8 :: sb, variance, detR, lh, sbO, sbmev
      real*8, parameter :: h = 1e-5, eps = 1e-13, eps2 = 1e-10
! eps avoid to become singular in 1st der & eps2 in 2nd der
      integer :: prev_ns, m_t, npx, counttimes

      contains

      Subroutine Setup_Kriging(nPoints,nInter,x_,dy_,y_)
#include "stdalloc.fh"
      integer :: nPoints, nInter, i, j
      real*8 :: x_(nInter,nPoints), dy_(nInter,nPoints), y_(nPoints)

      nInter_save = nInter
      nPoints_save = nPoints

      Call mma_Allocate(x,nInter,nPoints,Label="x")
      Call mma_Allocate(dy,nInter*nPoints,Label="dy")
      Call mma_Allocate(y,nPoints,Label="y")

!x is the n-dimensional internal coordinates
      x(:,:) = x_(:,:)
      ! write(6,*) 'x',x
!y is the energy
      y(:) = y_(:)
      ! write(6,*) 'y',y
!dy it's a vector of Grad-y (eq. (5)  ref. gradients of
! the energy with respect to the internal coordinates
      do i=1,nInter
        do j=1,nPoints
          dy((i-1)*nPoints+j) = dy_(i,j)
        enddo
      enddo

      return
      end subroutine Setup_kriging

      end module kriging_mod
