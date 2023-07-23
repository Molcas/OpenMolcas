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

use Definitions, only: wp, iwp

implicit none
private

public :: Diag, Diag_G, Diag_G_Hidden, Diag_Hidden, iiBstRSh, iiBstRSh_G, iiBstRSh_Hidden, iiBstRSh_L_Hidden, IndRed, IndRed_G, &
          IndRed_G_Hidden, IndRed_Hidden, IndRSh, IndRSh_G, IndRSh_G_Hidden, IndRSh_Hidden, InfRed, InfRed_G, InfRed_G_Hidden, &
          InfRed_Hidden, InfVec, InfVec_Bak, InfVec_G, InfVec_G_Hidden, InfVec_Hidden, iQuAB, iQuAB_here, iQuAB_Hidden, iQuAB_L, &
          iQuAB_L_Hidden, nnBstRSh, nnBstRSh_G, nnBstRSh_Hidden, nnBstRSh_L_Hidden, pTemp, pTemp1, pTemp3

integer(kind=iwp), allocatable, target :: iiBstRSh_Hidden(:,:,:), iiBstRSh_L_Hidden(:,:,:), IndRed_G_Hidden(:,:), &
                                          IndRed_Hidden(:,:), IndRSh_G_Hidden(:), IndRSh_Hidden(:), InfRed_G_Hidden(:), &
                                          InfRed_Hidden(:), InfVec_Bak(:,:,:), InfVec_G_Hidden(:,:,:), InfVec_Hidden(:,:,:), &
                                          iQuAB_here(:,:), iQuAB_Hidden(:,:), iQuAB_L_Hidden(:,:), nnBstRSh_Hidden(:,:,:), &
                                          nnBstRSh_L_Hidden(:,:,:)
integer(kind=iwp), pointer :: iiBstRSh(:,:,:) => null(), iiBstRSh_G(:,:,:) => null(), IndRed(:,:) => null(), &
                              IndRed_G(:,:) => null(), IndRSh(:) => null(), IndRSh_G(:) => null(), InfRed(:) => null(), &
                              InfRed_G(:) => null(), InfVec(:,:,:) => null(), InfVec_G(:,:,:) => null(), iQuAB(:,:) => null(), &
                              iQuAB_L(:,:) => null(), nnBstRSh(:,:,:) => null(), nnBstRSh_G(:,:,:) => null(), &
                              pTemp(:,:) => null(), pTemp1(:) => null(), pTemp3(:,:,:) => null()
real(kind=wp), allocatable, target :: Diag_G_Hidden(:), Diag_Hidden(:)
real(kind=wp), pointer :: Diag(:) => null(), Diag_G(:) => null()

end module ChoSwp
