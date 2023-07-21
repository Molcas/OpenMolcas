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

module ChoSwp

implicit none
private

public :: iQuAB, iQuAB_L, iQuAB_Hidden, iQuAB_L_Hidden, pTemp, iQuAB_here, nnBstRSh, nnBstRSh_Hidden, nnBstRSh_G, &
          nnBstRSh_L_Hidden, pTemp3, iiBstRSh, iiBstRSh_Hidden, iiBstRSh_G, iiBstRSh_L_Hidden, IndRSh, IndRSh_Hidden, IndRSh_G, &
          IndRSh_G_Hidden, pTemp1, InfRed, InfRed_Hidden, InfRed_G, InfRed_G_Hidden, InfVec, InfVec_Hidden, InfVec_G, &
          InfVec_G_Hidden, InfVec_Bak, IndRed, IndRed_Hidden, IndRed_G, IndRed_G_Hidden, Diag, Diag_Hidden, Diag_G, Diag_G_Hidden

integer, allocatable, target :: iQuAB_Hidden(:,:), iQuAB_L_Hidden(:,:), iQuAB_here(:,:)
integer, pointer :: iQuAB(:,:) => null(), iQuAB_L(:,:) => null(), pTemp(:,:) => null()

integer, allocatable, target :: nnBstRSh_Hidden(:,:,:), nnBstRSh_L_Hidden(:,:,:)
integer, allocatable, target :: iiBstRSh_Hidden(:,:,:), iiBstRSh_L_Hidden(:,:,:)
integer, pointer :: nnBstRSh(:,:,:) => null(), nnBstRSh_G(:,:,:) => null(), pTemp3(:,:,:) => null()
integer, pointer :: iiBstRSh(:,:,:) => null(), iiBstRSh_G(:,:,:) => null()
integer, pointer :: IndRSh(:) => null(), IndRSh_G(:) => null(), pTemp1(:) => null()
integer, allocatable, target :: IndRSh_Hidden(:), IndRSh_G_Hidden(:)
integer, pointer :: InfRed(:) => null(), InfRed_G(:) => null()
integer, allocatable, target :: InfRed_Hidden(:), InfRed_G_Hidden(:)
integer, pointer :: InfVec(:,:,:) => null(), InfVec_G(:,:,:) => null()
integer, allocatable, target :: InfVec_Hidden(:,:,:), InfVec_G_Hidden(:,:,:)
integer, allocatable, target :: InfVec_Bak(:,:,:)
integer, pointer :: IndRed(:,:) => null(), IndRed_G(:,:) => null()
integer, allocatable, target :: IndRed_Hidden(:,:), IndRed_G_Hidden(:,:)
real*8, pointer :: Diag(:) => null(), Diag_G(:) => null()
real*8, allocatable, target :: Diag_Hidden(:), Diag_G_Hidden(:)

end module ChoSwp
