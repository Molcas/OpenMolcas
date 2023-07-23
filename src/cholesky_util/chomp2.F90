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

use Definitions, only: wp, iwp

implicit none
private

type Pointer_2D
  real(kind=wp), pointer :: A(:,:) => null()
end type Pointer_2D

logical(kind=iwp) :: ChoMP2_allocated = .false., ChoMP2g_allocated = .false.
integer(kind=iwp), allocatable :: AdrR1(:,:,:), AdrR2(:,:,:), iFirst(:), iFirstS(:,:), LiMatij(:,:,:), LiPQprod(:,:,:), &
                                  LiT1am(:,:,:), LnBatOrb(:,:), LnMatij(:,:), LnOcc(:,:), LnPQprod(:,:), LnT1am(:,:), lUnit(:,:), &
                                  NumBatOrb(:), NumOcc(:)
real(kind=wp), allocatable, target :: MP2D_e_full(:), MP2D_full(:), MP2W_e_full(:), MP2W_full(:)
real(kind=wp), allocatable :: EFrozT(:), EOccuT(:), EVirtT(:), OldVec(:)
type(Pointer_2D) :: MP2D(8), MP2D_e(8), MP2W(8), MP2W_e(8)

public :: AdrR1, AdrR2, ChoMP2_allocated, ChoMP2g_allocated, EFrozT, EOccuT, EVirtT, iFirst, iFirstS, LiMatij, LiPQprod, LiT1am, &
          LnBatOrb, LnMatij, LnOcc, LnPQprod, LnT1am, lUnit, MP2D, MP2D_e, MP2D_e_full, MP2D_full, MP2W, MP2W_e, MP2W_e_full, &
          MP2W_full, NumBatOrb, NumOcc, OldVec, Pointer_2D

end module ChoMP2
