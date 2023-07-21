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

module ChoMP2

implicit none
private

public :: Pointer_2D
public :: ChoMP2_allocated, ChoMP2g_allocated, OldVec, EFrozT, EOccuT, EVirtT
public :: AdrR1, AdrR2
public :: MP2W_full, MP2W
public :: MP2D_full, MP2D
public :: MP2W_e_full, MP2W_e
public :: MP2D_e_full, MP2D_e
public :: iFirst, iFirstS, NumOcc, LnOcc, LnT1am, LiT1am, LnMatij, LiMatij
public :: lUnit, NumBatOrb, LnBatOrb
public :: LnPQprod, LiPQprod

logical :: ChoMP2_allocated = .false.
logical :: ChoMP2g_allocated = .false.
real*8, allocatable :: OldVec(:)

type Pointer_2D
  real*8, pointer :: A(:,:) => null()
end type Pointer_2D

real*8, allocatable, target :: MP2D_full(:)
type(Pointer_2D) :: MP2D(8)
real*8, allocatable, target :: MP2W_full(:)
type(Pointer_2D) :: MP2W(8)
real*8, allocatable, target :: MP2D_e_full(:)
type(Pointer_2D) :: MP2D_e(8)
real*8, allocatable, target :: MP2W_e_full(:)
type(Pointer_2D) :: MP2W_e(8)

real*8, allocatable :: EFrozT(:)
real*8, allocatable :: EOccuT(:)
real*8, allocatable :: EVirtT(:)

integer, allocatable :: AdrR1(:,:,:), AdrR2(:,:,:)

integer, allocatable :: iFirst(:), iFirstS(:,:), NumOcc(:), LnOcc(:,:)
integer, allocatable :: LnT1am(:,:), LiT1am(:,:,:)
integer, allocatable :: LnMatij(:,:), LiMatij(:,:,:)
integer, allocatable :: lUnit(:,:)
integer, allocatable :: NumBatOrb(:), LnBatOrb(:,:)
integer, allocatable :: LnPQprod(:,:), LiPQprod(:,:,:)

end module ChoMP2
