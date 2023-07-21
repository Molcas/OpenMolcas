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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************

module ChoArr

implicit none
private

integer, allocatable :: iSOShl(:)
integer, allocatable :: iShlSO(:)
integer, allocatable :: iBasSh(:,:)
integer, allocatable :: nBasSh(:,:)
integer, allocatable :: nBstSh(:)
integer, allocatable :: iSP2F(:)
integer, allocatable :: iAtomShl(:)
integer, allocatable :: iRS2F(:,:)
integer, allocatable :: IntMap(:)
integer, allocatable :: iScr(:)
integer, allocatable :: nDimRS(:,:)
integer, allocatable :: iL2G(:)
integer, allocatable :: iShP2RS(:,:)
integer, allocatable :: iShP2Q(:,:)
integer, allocatable :: iSimRI(:)

integer, allocatable :: iOff_Batch(:,:)
integer :: nDim_Batch(8)

integer, allocatable :: iQL2G(:,:)

type rPointers
  real*8, pointer :: Array(:,:) => null()
end type rPointers

type(rPointers) :: LQ(8)
real*8, allocatable, target :: LQ_Tot(:)

integer :: nQual_L(8)

integer, allocatable :: Idle(:)
integer, allocatable :: MySP(:)
integer :: n_MySP

public :: iSOShl, iShlSO, iBasSh, nBasSh, nBstSh, iSP2F, iAtomShl, iRS2F, IntMap, iScr, nDimRS, iL2G, iShP2RS, iShP2Q, iOff_Batch, &
          nDim_Batch, iQL2G, LQ, LQ_Tot, nQual_L, Idle, MySP, n_MySP, iSimRI

end module ChoArr
