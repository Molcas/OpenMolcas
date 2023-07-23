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

use Definitions, only: wp, iwp

implicit none
private

type rPointers
  real(kind=wp), pointer :: Array(:,:) => null()
end type rPointers

integer(kind=iwp) :: n_MySP, nDim_Batch(8), nQual_L(8)
integer(kind=iwp), allocatable :: iAtomShl(:), iBasSh(:,:), Idle(:), iL2G(:), IntMap(:), iOff_Batch(:,:), iQL2G(:,:), iRS2F(:,:), &
                                  iScr(:), iShlSO(:), iShP2Q(:,:), iShP2RS(:,:), iSimRI(:), iSOShl(:), iSP2F(:), MySP(:), &
                                  nBasSh(:,:), nBstSh(:), nDimRS(:,:)
real(kind=wp), allocatable, target :: LQ_Tot(:)
type(rPointers) :: LQ(8)

public :: iAtomShl, iBasSh, Idle, iL2G, IntMap, iOff_Batch, iQL2G, iRS2F, iScr, iShlSO, iShP2Q, iShP2RS, iSimRI, iSOShl, iSP2F, &
          LQ, LQ_Tot, MySP, n_MySP, nBasSh, nBstSh, nDim_Batch, nDimRS, nQual_L

end module ChoArr
