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
! Copyright (C) 2000, Roland Lindh                                     *
!***********************************************************************

function Compute_Tau(Weights,mGrid,iSpin)
!***********************************************************************
!      Author:Roland Lindh, Department of Chemical Physics, University *
!             of Lund, SWEDEN. November 2000                           *
!***********************************************************************

use nq_Grid, only: Tau
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Compute_Tau
integer(kind=iwp), intent(in) :: mGrid, iSpin
real(kind=wp), intent(in) :: Weights(mGrid)
integer(kind=iwp) :: iGrid
real(kind=wp) :: TauA

!                                                                      *
!***********************************************************************
!                                                                      *
Compute_Tau = Zero
if (iSpin == 1) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! iSpin == 1

  do iGrid=1,mGrid

    TauA = Two*Tau(1,iGrid)

    ! Accumulate contributions to the integrated Tau

    Compute_Tau = Compute_Tau+TauA*Weights(iGrid)

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! iSpin /= 1

  do iGrid=1,mGrid

    TauA = (Tau(1,iGrid)+Tau(2,iGrid))

    ! Accumulate contributions to the integrated density

    Compute_Tau = Compute_Tau+TauA*Weights(iGrid)

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *

return

end function Compute_Tau
